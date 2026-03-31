# Uruk-Hai GPU Collision Review — gpu_collision.cpp

## Summary
Tested against: GL resource leaks, bit-packing alignment, null pointer dereference, shader binding order, glGetError coverage.

**CONFIRMED BUGS: 3**

---

## BUG #1: GL Resource Leak on init_context() Failure

**Location**: `GpuCollisionEvaluator::init_context()` lines 249–327

**Severity**: CRITICAL resource leak

**Issue**: If `wglCreateContext()` fails at line 293, the code returns false but does not call `ReleaseDC()` on the device context acquired at line 273. The DC remains leaked.

```cpp
HDC hdc = GetDC(hwnd);
if (!hdc) {
    BOOST_LOG_TRIVIAL(warning) << "Snuggle GPU: GetDC failed";
    return false;
}
gl_dc_ = (void*)hdc;

// ... later ...

HGLRC ctx = wglCreateContext(hdc);
if (!ctx) {
    BOOST_LOG_TRIVIAL(warning) << "Snuggle GPU: wglCreateContext failed";
    return false;  // ← DC NOT RELEASED
}
```

Also fails to release DC if `glewInit()` fails (line 308) and if OpenGL version check fails (line 321).

**Fix**: Store hwnd and release on all error paths, or wrap in RAII guard. At minimum:
```cpp
if (!ctx) {
    ReleaseDC(hwnd, hdc);
    return false;
}
```

---

## BUG #2: HDC Not Released on Multiple Error Paths

**Location**: `GpuCollisionEvaluator::init_context()` lines 293–332

**Severity**: CRITICAL resource leak

**Issue**: In addition to Bug #1, the DC is also not released when:
- `glewInit()` fails (line 308–311)
- OpenGL version is too old (line 321–325)

Both error paths return without calling `ReleaseDC(hwnd, hdc)`.

**Fix**: Consolidate cleanup. Use a local variable to track the dc:
```cpp
HDC hdc = GetDC(hwnd);
if (!hdc) {
    BOOST_LOG_TRIVIAL(warning) << "Snuggle GPU: GetDC failed";
    DestroyWindow(hwnd);
    return false;
}

// ... all subsequent errors must:
ReleaseDC(hwnd, hdc);
DestroyWindow(hwnd);  // if hwnd still needed
return false;
```

---

## BUB #3: Missing glGetError Check After glShaderSource + glCompileShader

**Location**: `GpuCollisionEvaluator::compile_shader()` lines 335–374

**Severity**: MEDIUM (logic still works but errors hidden)

**Issue**: After `glCompileShader()` at line 344, there is no `glGetError()` check. If the shader source pointer is invalid (nullptr, despite the const string), the compile silently fails in GL but returns GL_TRUE on `glGetShaderiv(...GL_COMPILE_STATUS...)` in some drivers. The code relies on the compile status flag alone, which is correct, but does not detect GL state corruption.

This is actually safe-by-accident here because glGetShaderiv() is called immediately after, and status is checked. However, best practice would check for GL errors after shader operations:

```cpp
glShaderSource(shader, 1, &COLLISION_SHADER_SRC, nullptr);
glCompileShader(shader);
GLenum err = glGetError();  // <- Missing
if (err != GL_NO_ERROR) {
    BOOST_LOG_TRIVIAL(warning) << "Snuggle GPU: GL error after compile: 0x" << std::hex << err;
    // still check compile status below, but now we know about GL state errors
}
```

---

## BUG #4: Null Dereference Risk in cleanup() on Platform Mismatch

**Location**: `GpuCollisionEvaluator::cleanup()` lines 549–580

**Severity**: LOW (Windows-only code path, unlikely in practice)

**Issue**: `cleanup()` unconditionally calls `wglMakeCurrent(hdc, ctx)` at line 556 if `gl_context_` is set, but does not verify that `gl_dc_` is non-null. On non-Windows platforms where `init_context()` returns early (line 330), `gl_dc_` remains nullptr. If someone compiles for a hypothetical platform and `gl_context_` gets set somehow, calling `wglMakeCurrent(NULL, ...)` could have platform-specific behavior.

```cpp
void GpuCollisionEvaluator::cleanup()
{
    if (gl_context_) {
#ifdef _WIN32
        HDC hdc = (HDC)gl_dc_;  // ← gl_dc_ could be nullptr if init_context() failed early
        HGLRC ctx = (HGLRC)gl_context_;
        wglMakeCurrent(hdc, ctx);  // ← Passes nullptr hdc
#endif
```

This is mitigated by the fact that `init_context()` only sets `gl_context_` after successfully getting `gl_dc_` (line 298), but the logic is fragile.

**Fix**: Add defensive check:
```cpp
if (gl_context_ && gl_dc_) {
```

---

## Non-Bugs (Verified Safe)

### Bit-Packing Alignment (lines 405–437): CLEAN
The uint8→uint32 repacking is correct:
- Source: bit N at `src[N/8] >> (N%8)`
- Dest: bit N at `dest_words[N/32] |= (1u << (N%32))`
- Both access the linear bit array in the same order. The destination is zero-initialized (line 402), so |= accumulation is safe.
- Byte alignment: `offset` is incremented by `grid_bytes = ((total_bits+31)/32)*4`, which maintains 4-byte alignment on each grid boundary.

### Shader Uniform Setting (lines 514–518): CLEAN
Uniforms are set AFTER `glUseProgram(program_)` at line 514. The order is correct. All glUniform* calls (515–518) reference `program_`, which is now current.

### glGetError Checks: MOSTLY GOOD
- Lines 461–468: ✓ Checked after upload_grids
- Lines 540–546: ✓ Checked after evaluate_batch
- Missing after compile_shader (noted as Bug #3 above, but safe)

### SSBO Binding (lines 508–511): CLEAN
All four glBindBufferBase calls use sequential bindings 0–3, matching shader layout declarations (lines 131–134).

### Results Readback (lines 525–535): CLEAN
- glGetBufferSubData implicitly syncs (GPU-CPU), so no explicit barrier needed after glMemoryBarrier.
- Individual struct packing is correct (16 bytes for PlacementGPU, verified by static_assert at line 109).
- Grid metadata packing is correct (48 bytes, verified by static_assert at line 102).

---

## Recommendations

1. **Immediate** (Bug #1, #2): Add DC release on all error paths in init_context().
2. **Defensive** (Bug #4): Add null check for gl_dc_ in cleanup().
3. **Nice-to-have** (Bug #3): Add glGetError check after glCompileShader for consistency.

