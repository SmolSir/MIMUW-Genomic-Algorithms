import numpy as np

# Initialize matrix D with given dimensions and fill the first row and first column
rows, cols = 12, 20
edist = 3
D = np.zeros((rows + 1, cols + 1)).astype(int)
# D[1:, 0] = np.arange(1, rows + 1)
# D[1:, cols] = np.arange(cols + 1, cols + rows + 1)
# D[0, 1:] = np.arange(1, cols + 1)
# D[rows, 1:] = np.arange(rows + 1, rows + cols + 1)

bruh = np.random.randint(0, 2, (rows + 1, cols + 1))

print(f"matrix: \n{D}")
# Compute each diagonal from top-left to bottom-right


left_bottom_bound_diag_row_start = edist
left_bottom_bound_diag_row_end = rows + 1

left_bottom_bound_diag_col_start = 0
left_bottom_bound_diag_col_end = max(0, rows - edist + 1)

left_bottom_bound_diag_row_indices = np.arange(left_bottom_bound_diag_row_start, left_bottom_bound_diag_row_end)
left_bottom_bound_diag_col_indices = np.arange(left_bottom_bound_diag_col_start, left_bottom_bound_diag_col_end)

D[left_bottom_bound_diag_row_indices, left_bottom_bound_diag_col_indices] = edist * 10

top_right_bound_diag_row_start = 0
top_right_bound_diag_row_end = max(0, rows - edist + 1)

top_right_bound_diag_col_start = max(0, cols - rows + edist)
top_right_bound_diag_col_end = cols + 1

top_right_bound_diag_row_indices = np.arange(top_right_bound_diag_row_start, top_right_bound_diag_row_end)
top_right_bound_diag_col_indices = np.arange(top_right_bound_diag_col_start, top_right_bound_diag_col_end)

D[top_right_bound_diag_row_indices, top_right_bound_diag_col_indices] = edist * 10

print(f"matrix with blockers: \n{D}")


for diag in range(2, rows + cols + 1):
    print(f"diag: {diag}")
    # Determine the start and end for row and column indices along the diagonal

    row_start = min(rows, diag - 1 if diag <= edist else edist + (diag - edist - 1) // 2)
    row_end = max(0, (diag - edist) // 2 - 4 if diag <= cols + rows - edist else diag - cols - 1)

    col_start = max(1, diag - rows if diag > 2 * rows - edist + 2 else (diag - edist) // 2 + 1)
    col_end = min(cols + 1, diag if diag <= cols - rows + edist else cols - rows + edist + (diag - cols + rows - edist + 1) // 2)

    # Get the indices for this diagonal
    row_indices = np.arange(row_start, row_end, -1)
    col_indices = np.arange(col_start, col_end)

    # Filter indices to stay within the bounds of the matrix
    # row_indices = row_indices[row_indices < diag]
    # col_indices = col_indices[col_indices < cols]

    print(f"row_indices: {row_indices}")
    print(f"col_indices: {col_indices}")
    print(f"diag_values: {D[row_indices, col_indices]}")
    # Calculate the minimum for each cell on this diagonal
    # D[row_indices, col_indices] = diag
    D[row_indices, col_indices] = np.minimum.reduce([
        D[row_indices - 1, col_indices - 1] + bruh[row_indices - 1, col_indices - 1], # Top-left
        D[row_indices - 1, col_indices] + 1,    # Top
        D[row_indices, col_indices - 1] + 1    # Left
    ])

print("\nComputed matrix:")
print(D)
