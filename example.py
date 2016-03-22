from bttt.model import import_tree_model
from bttt.trees import get_ts, DeterministicTree





tree = DeterministicTree(40)

for n in tree.nodes:
    tree.values[n] = 0.1


model = import_tree_model('example.yaml', tree=tree)

model

sol = model.solve()


sim = get_ts(tree, sol, 'e')


print(sim)
