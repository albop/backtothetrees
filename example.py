from bttt.model import import_tree_model
from bttt.trees import get_ts, DeterministicTree



etree = DeterministicTree(N)
for s in etree.nodes:
    etree.values[s] = shock

for n in tree.nodes:
    tree.values[n] = 0.1


model = import_tree_model('example.yaml', tree=tree)

model

sol = model.solve()


sim = get_ts(tree, sol, 'e')


print(sim)
