import os

# Directory containing the files
directory = r'C:\GitHub\kwave\Matlab_codes\newdatafinal\result2500'

# Total number of expected files
total_files = 2500

# List to hold missing file indices
missing_indices = []

# Check each expected file
for i in range(1, total_files + 1):
    file_name = 'dataset_{}.mat'.format(i)
    if not os.path.isfile(os.path.join(directory, file_name)):
        missing_indices.append(i)

# Save missing indices to a file
with open('missing_files.txt', 'w') as f:
    for index in missing_indices:
        f.write('{}\n'.format(index))

print('Missing files list saved to missing_files.txt')
