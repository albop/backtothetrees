from bttt.model import import_tree_model
from bttt.trees import get_ts, DeterministicTree

N = 40
etree = DeterministicTree(N)
for n in etree.nodes:
    etree.values[n] = 0.1

model = import_tree_model('example.yaml', tree=etree)

sol = model.solve()

sim = get_ts(etree, sol, 'e')


print(sim)
