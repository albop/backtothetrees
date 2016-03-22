A model begins with the definition of an event tree, that is a set of nodes
linked by probabilities of transitions, without going back.

```
from bttt import DeterministicTree
tree = DeterministicTree(10) # 10 periods
```

We set the exogenous values at each node of the tree:

```python
for n in tree.nodes:
    tree.values[n] = 0.1
```

A model, on a tree can be defined a yaml file. For instance `example.yaml` contains a section:

```yaml
equations:
  - -z[s] + etree.values[s]                                         | z[s]
  - -e[s] + a/(a+c)*E[ e[x]   | x in S(s)] + 1.0/(a+c)*(z[s]-f[s])  | e[s]
  - len(s)==1:
    - -Gamma[s] + (e[s]-estar)                                      | Gamma[s]
  - len(s)>1:
    - -Gamma[s] + a/(a+c)*Gamma[P(s)] + beta**t*(e[s]-estar)        | Gamma[s]
  - len(s)==1:
    - -Gamma[s] + E[ Gamma[x] | x in S(s) ] + psi[s] - phi[s]       | f[s]
  - 1<len(s)<len(etree):
    - -Gamma[s] + E[ Gamma[x] | x in S(s) ] + psi[s] - phi[s]       | f[s]
  - len(s)==len(etree):
    - -f[s]                                                         | f[s]
  - -psi[s] + 1000*( Sum[f[x] | x in H(s)] - Rbar  )                | 0 <= psi[s]
  - -phi[s] + (-(f[s]-min_f))*1000                                  | 0 <= phi[s]
```

By convention, when `z` is a variable, `z[s]` denotes the value of `z` at a node
`s` of the event tree. Equations can be applied to a subset of the nodes, by preceding
the equation by a condition, for instance `1<len(s)<len(etree)`.
Complementarity conditions can optionnally be associated to equations.
Remark the following elements of syntax:
- `P(s)` denotes the predecessor of node `s`
- `H(s)` denotes the set of all ancestors of node `s`
- `S(s)` denotes the set of all successors of node `s` with the probabilities of transitions.
- `Sum[expr[x] | x in H(s)]` computes the sum of the `expr[x]` where `x` exhausts the set of all ancestors of `s`.
- `E[ Gamma[x] | x in S(s) ]` is the expectation of `Gamma[s]` over all successors of `s`


The calibration of the model is provided in another section:

```yaml
calibration:
    beta: 0.8421052631578947
    a: 0.8
    c: 0.15
    estar: -0.0
    Rbar: 1.0
    min_f: 0
    kappa: 1.0
    N: 40
    zbar: 0.1
    p: 1
```

It is then possible to import and solve the model with:

```python
model = import_tree_model('example.yaml', tree=dt)
sol = model.solve()
```

The result `sol` is a mapping associating all variables at all nodes with a value.
It can be used to compute simulations.

```python
sim = get_ts(tree, sol, 'e')
print(sim)
```
