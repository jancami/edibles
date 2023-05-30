# sightlines = ['abc', 'cde', 'efg']
# B = 0.01
# T = 75.5
# delta_B = -0.1
# zeta = -0.1 
# sigma = 0.02
# origin =  0

# T_values = [None] * len(sightlines) 
# sigma_values = [None] * len(sightlines) 
# origin_values = [None] * len(sightlines) 

# my_list = [B, delta_B, zeta] + T_values + sigma_values + origin_values

# print(my_list)

# first_T_index = 3
# last_T_index = first_T_index + len(T_values) 

# first_sigma_index = last_T_index  
# last_sigma_index = first_sigma_index + len(sigma_values)

# first_origin_index = last_sigma_index  
# last_origin_index = first_origin_index + len(origin_values)

# #print(first_sigma_index)

# for t in range(first_T_index, last_T_index ):
#     print(t)
#     my_list[t] = T

# for s in range(first_sigma_index, last_sigma_index):
#     print(s)
#     my_list[s] = sigma
    
# for o in range(first_origin_index, last_origin_index):
#     print(o)
#     my_list[o] = origin
    
    
# print(my_list)

sightlines = ['abc', 'cde']
B = 0.01
T = 75.5
delta_B = -0.1
zeta = -0.1 
sigma = 0.02
origin =  0

sightlines = ['166937', '000000']

params_list = [B, delta_B, zeta]
for i in range(len(sightlines)):
    exec(f"T{i+1} = T")
    params_list.append(eval(f"T{i+1}"))
   
    
for i in range(len(sightlines)):    
    exec(f"sigma{i+1} = sigma")    
    params_list.append(eval(f"sigma{i+1}"))
    
for i in range(len(sightlines)):    
    exec(f"origin{i+1} = origin")    
    params_list.append(eval(f"origin{i+1}"))

print(params_list)

first_T_index = 3
last_T_index = first_T_index + len(sightlines) 

first_sigma_index = last_T_index  
last_sigma_index = first_sigma_index + len(sightlines) 

first_origin_index = last_sigma_index  
last_origin_index = first_origin_index +len(sightlines) 

T_values = params_list[first_T_index:last_T_index]
sigma_values = params_list[first_sigma_index:last_sigma_index]
origin_values = params_list[first_origin_index:last_origin_index]

# print(T_values)
# print(sigma_values)
# print(origin_values)

# wave_list = [1,2,3,5,6,7]
# print(wave_list[-1])

# list1 = [-2, -1, 0, 1, 3]
# list2 = [-1, 0, 1, 2]

# increment = list1[-1] + 1  # Determine the increment value

# new_list = list1 + [x + increment for x in list2]
# print(new_list)

