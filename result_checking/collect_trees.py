
part_file = 'wu_partitionout.nex'

part_list = []
with open(part_file) as partition:
        for line in partition.readlines():
            if '=' in line:
                list_line = line.split()
                part_list.append(list_line[0])

output_file = 'trees_1c.txt'
n = 0
with open(output_file, 'w') as output_file:
    for i in range(len(part_list)):
        tree_file = 'No' + str(i+1) + part_list[i] + '_r/logs/1class.treefile'
        with open(tree_file, 'r') as current_file:
            output_file.write(current_file.read())
            n = n + 1
print(n)


output_file = 'trees_2c.txt'
n = 0
with open(output_file, 'w') as output_file:
    for i in range(len(part_list)):
        tree_file = 'No' + str(i+1) + part_list[i] + '_r/logs/fix_model.treefile'
        with open(tree_file, 'r') as current_file:
            output_file.write(current_file.read())
            n = n + 1
print(n)
