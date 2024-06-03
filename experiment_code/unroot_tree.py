from ete3 import Tree

true_treefile = r'C:\Users\u7151703\Desktop\research\results\simulation\mast_evo\astral\cannon_1c.tre'
true_treestr = open(true_treefile,'r').read()

true_tree = Tree(true_treestr,format = 1)
true_tree.unroot()
with open(r'C:\Users\u7151703\Desktop\research\results\simulation\mast_evo\astral\cannon_1c_unrooted.treefile', 'w+') as result:
    result.write(true_tree.write(format = 5) + '\n')    


true_treefile = r'C:\Users\u7151703\Desktop\research\results\simulation\mast_evo\astral\cannon_2c.tre'
true_treestr = open(true_treefile,'r').read()

true_tree = Tree(true_treestr,format = 1)
true_tree.unroot()
with open(r'C:\Users\u7151703\Desktop\research\results\simulation\mast_evo\astral\cannon_2c_unrooted.treefile', 'w+') as result:
    result.write(true_tree.write(format = 5) + '\n')    


true_treefile = r'C:\Users\u7151703\Desktop\research\results\simulation\mast_evo\astral\ran_1c.tre'
true_treestr = open(true_treefile,'r').read()

true_tree = Tree(true_treestr,format = 1)
true_tree.unroot()
with open(r'C:\Users\u7151703\Desktop\research\results\simulation\mast_evo\astral\ran_1c_unrooted.treefile', 'w+') as result:
    result.write(true_tree.write(format = 5) + '\n')    


true_treefile = r'C:\Users\u7151703\Desktop\research\results\simulation\mast_evo\astral\ran_2c.tre'
true_treestr = open(true_treefile,'r').read()

true_tree = Tree(true_treestr,format = 1)
true_tree.unroot()
with open(r'C:\Users\u7151703\Desktop\research\results\simulation\mast_evo\astral\ran_2c_unrooted.treefile', 'w+') as result:
    result.write(true_tree.write(format = 5) + '\n')    


true_treefile = r'C:\Users\u7151703\Desktop\research\results\simulation\mast_evo\astral\wu_1c.tre'
true_treestr = open(true_treefile,'r').read()

true_tree = Tree(true_treestr,format = 1)
true_tree.unroot()
with open(r'C:\Users\u7151703\Desktop\research\results\simulation\mast_evo\astral\wu_1c_unrooted.treefile', 'w+') as result:
    result.write(true_tree.write(format = 5) + '\n')    


true_treefile = r'C:\Users\u7151703\Desktop\research\results\simulation\mast_evo\astral\wu_2c.tre'
true_treestr = open(true_treefile,'r').read()

true_tree = Tree(true_treestr,format = 1)
true_tree.unroot()
with open(r'C:\Users\u7151703\Desktop\research\results\simulation\mast_evo\astral\wu_2c_unrooted.treefile', 'w+') as result:
    result.write(true_tree.write(format = 5) + '\n')    



