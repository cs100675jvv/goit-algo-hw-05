def binary_search(arr, target):
    left, right = 0, len(arr) - 1
    iterations = 0
    upper_bound = None

    while left <= right:
        iterations += 1
        mid = (left + right) // 2

        if arr[mid] == target:
            return (iterations, arr[mid])
        elif arr[mid] < target:
            left = mid + 1
        else:
            upper_bound = arr[mid]
            right = mid - 1

    return (iterations, upper_bound)

# Приклад використання
sorted_array = [0.1, 0.3, 0.5, 1.1, 1.5, 2.0, 2.3, 2.9, 3.4]
target_value = 1.1

result = binary_search(sorted_array, target_value)
print(result)  
