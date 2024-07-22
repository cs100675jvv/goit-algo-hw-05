import timeit
# import re

# Read the content of the text files
with open('task_3\стаття 1.txt', 'r', encoding='utf-8') as file:
    text1 = file.read()

with open('task_3/стаття 2.txt', 'r', encoding='utf-8') as file:
    text2 = file.read()

# Substring to search (one existing, one non-existing)
existing_substring = "алгоритми"
non_existing_substring = "неіснуючийпідрядок"

# Boyer-Moore Algorithm Implementation
def boyer_moore(text, pattern):
    m = len(pattern)
    n = len(text)
    if m == 0:
        return 0
    last = {}
    for i in range(m):
        last[pattern[i]] = i
    i = m - 1
    k = m - 1
    while i < n:
        if text[i] == pattern[k]:
            if k == 0:
                return i
            else:
                i -= 1
                k -= 1
        else:
            j = last.get(text[i], -1)
            i = i + m - min(k, j + 1)
            k = m - 1
    return -1

# Knuth-Morris-Pratt Algorithm Implementation
def kmp_search(text, pattern):
    n = len(text)
    m = len(pattern)
    lps = [0] * m
    j = 0
    compute_lps_array(pattern, m, lps)
    i = 0
    while i < n:
        if pattern[j] == text[i]:
            i += 1
            j += 1
        if j == m:
            return i - j
        elif i < n and pattern[j] != text[i]:
            if j != 0:
                j = lps[j-1]
            else:
                i += 1
    return -1

def compute_lps_array(pattern, m, lps):
    length = 0
    lps[0] = 0
    i = 1
    while i < m:
        if pattern[i] == pattern[length]:
            length += 1
            lps[i] = length
            i += 1
        else:
            if length != 0:
                length = lps[length-1]
            else:
                lps[i] = 0
                i += 1

# Rabin-Karp Algorithm Implementation
def rabin_karp(text, pattern):
    d = 256
    q = 101
    m = len(pattern)
    n = len(text)
    p = 0
    t = 0
    h = 1
    for i in range(m-1):
        h = (h * d) % q
    for i in range(m):
        p = (d * p + ord(pattern[i])) % q
        t = (d * t + ord(text[i])) % q
    for i in range(n - m + 1):
        if p == t:
            for j in range(m):
                if text[i + j] != pattern[j]:
                    break
            j += 1
            if j == m:
                return i
        if i < n - m:
            t = (d * (t - ord(text[i]) * h) + ord(text[i + m])) % q
            if t < 0:
                t = t + q
    return -1

# Measure execution time for each algorithm with the existing substring
time_existing_bm = timeit.timeit(lambda: boyer_moore(text1, existing_substring), number=1000)
time_existing_kmp = timeit.timeit(lambda: kmp_search(text1, existing_substring), number=1000)
time_existing_rk = timeit.timeit(lambda: rabin_karp(text1, existing_substring), number=1000)

time_non_existing_bm = timeit.timeit(lambda: boyer_moore(text1, non_existing_substring), number=1000)
time_non_existing_kmp = timeit.timeit(lambda: kmp_search(text1, non_existing_substring), number=1000)
time_non_existing_rk = timeit.timeit(lambda: rabin_karp(text1, non_existing_substring), number=1000)

# Print results for text1
print("Text 1 - Existing Substring:")
print(f"Boyer-Moore: {time_existing_bm}")
print(f"KMP: {time_existing_kmp}")
print(f"Rabin-Karp: {time_existing_rk}")

print("Text 1 - Non-Existing Substring:")
print(f"Boyer-Moore: {time_non_existing_bm}")
print(f"KMP: {time_non_existing_kmp}")
print(f"Rabin-Karp: {time_non_existing_rk}")

# Measure execution time for each algorithm with the existing substring in text2
time_existing_bm_2 = timeit.timeit(lambda: boyer_moore(text2, existing_substring), number=1000)
time_existing_kmp_2 = timeit.timeit(lambda: kmp_search(text2, existing_substring), number=1000)
time_existing_rk_2 = timeit.timeit(lambda: rabin_karp(text2, existing_substring), number=1000)

time_non_existing_bm_2 = timeit.timeit(lambda: boyer_moore(text2, non_existing_substring), number=1000)
time_non_existing_kmp_2 = timeit.timeit(lambda: kmp_search(text2, non_existing_substring), number=1000)
time_non_existing_rk_2 = timeit.timeit(lambda: rabin_karp(text2, non_existing_substring), number=1000)

# Print results for text2
print("Text 2 - Existing Substring:")
print(f"Boyer-Moore: {time_existing_bm_2}")
print(f"KMP: {time_existing_kmp_2}")
print(f"Rabin-Karp: {time_existing_rk_2}")

print("Text 2 - Non-Existing Substring:")
print(f"Boyer-Moore: {time_non_existing_bm_2}")
print(f"KMP: {time_non_existing_kmp_2}")
print(f"Rabin-Karp: {time_non_existing_rk_2}")
