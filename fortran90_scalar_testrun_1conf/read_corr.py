import sys
import numpy as np

filename = sys.argv[1]

with open(filename, "rb") as f:
    dims = np.fromfile(f, np.int32, 4)
    data = np.fromfile(f, np.complex128)

arr = data.reshape(tuple(dims), order="F")

print("shape =", dims.tolist())
print("first 20 values:")
print(arr.reshape(-1, order="F")[:20])
