
def draw_graph(graph, pos=None):
    from numpy import array
    import networkx as nx
    import bqplot as bq

    if pos is None:
        pos = nx.spring_layout(graph)
    nodes = graph.nodes()
    x = [pos[n][0] for n in nodes]
    y = [pos[n][1] for n in nodes]
    lines_x = []
    lines_y = []
    for k,v in graph.edges():
        i1 = graph.nodes().index(k)
        i2 = graph.nodes().index(v)
        lines_x.append([x[i1], x[i2]])
        lines_y.append([y[i1], y[i2]])
    lines_x = array(lines_x)
    lines_y = array(lines_y)

    x_sc = bq.LinearScale()
    y_sc = bq.LinearScale()
    points = bq.Scatter(x=x, y=y, scales={'x': x_sc, 'y': y_sc})
    lines = bq.Lines(x=lines_x,y=lines_y, scales={'x': x_sc, 'y': y_sc}, colors=['red'])
    ax_x = bq.Axis(scale=x_sc)
    ax_y = bq.Axis(scale=y_sc, orientation='vertical')
    fig = bq.Figure(marks=[lines,points], axes=[ax_x, ax_y])
    return fig


def draw_tree(tree,pos=None):
    if pos is None:
        pos = layout_tree(tree)
    return draw_graph(tree.graph(),pos=pos)

def layout_tree(tree):
    pos = {}
    for n in tree.nodes:
        x = len(n)
        y = -sum(n[:x])
        pos[n] = x,y
    return pos
