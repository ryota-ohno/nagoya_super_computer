array_1=[[0,0,0],[1,1,1]]
array_2=[[2,2,2],[3,3,3]]
array_list=[array_1,array_2]
concatenated = []
for arr in array_list:
    concatenated.extend(arr)
print(concatenated)