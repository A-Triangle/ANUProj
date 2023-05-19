import numpy as np
import json

with open('Trial.json', 'rt') as f:
        data = json.load(f)[0]

matrix = np.array(data['predicted_aligned_error'], dtype=np.float64)

i=0
while i <= matrix.shape[1]:
        print(matrix.argmin(axis=1))
        i+=1
        