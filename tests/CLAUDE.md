# Testing Guide for OrcaSlicer

Catch2 v2 (`#include <catch2/catch.hpp>`). For full reference (patterns, CMake, known issues): read `TESTING_REFERENCE.md`.

## Critical Rules

### 1. SECTIONS IN LOOPS — NEVER REUSE NAMES
```cpp
// WRONG — unpredictable behavior
SECTION("Same name") { ... }  // inside a loop

// CORRECT — unique per iteration
DYNAMIC_SECTION("Section " << i) { ... }
```

### 2. THREAD SAFETY — ASSERTIONS ARE NOT THREAD-SAFE
```cpp
// WRONG — undefined behavior
threads.emplace_back([]() { REQUIRE(x == y); });

// CORRECT — collect results, assert on main thread
std::atomic<int> passed{0};
threads.emplace_back([&]() { if (x == y) passed++; });
for (auto& t : threads) t.join();
REQUIRE(passed == expected_count);
```

### 3. EXPRESSION DECOMPOSITION — SPLIT BINARY OPERATORS
```cpp
// WRONG — shows "false" on failure, not values
REQUIRE(a > 0 && b < 10);

// CORRECT — each shows individual values
REQUIRE(a > 0);
REQUIRE(b < 10);
```

### 4. FLOATING POINT — USE MATCHERS, NOT APPROX
```cpp
// WRONG — Approx is deprecated
REQUIRE(val == Catch::Approx(expected));

// CORRECT
REQUIRE_THAT(val, WithinAbs(expected, 0.001));
REQUIRE_THAT(val, WithinRel(expected, 0.01));
```

### 5. TEST ORDERING — ALWAYS RANDOM
```bash
./tests --order rand --warn NoAssertions
```
