from dolang.codegen import to_source


import copy
import yaml
import time

from .pattern import ReplaceExpectation, match
from .symbols import LinearExpression
from .trees import *
from .trees import get_ts
from .symbols import LIVar

def import_tree_model(filename, tree=None, key='equations'):

    with open(filename) as f:
        d = yaml.load(f)

    lines = []
    for k in d[key]:
        if isinstance(k,str):
            lines.append(k)
        else:
            cond = next(iter(k.keys()))
            val = k[cond][0]
            s = cond + ' : ' + val
            lines.append(s)

    calibration = d['calibration']

    if tree is None:
        return lines

    else:
        return Model(tree, lines, calibration)

def read_equations(lines):

    lines = [l.strip() for l in lines if len(l.strip())>=1]
    conditions = []
    equations = []
    variables = []
    complementarities = []

    for l in lines:
        if ':' in l:
            cond, rhs = str.split(l, ':')
        else:
            cond, rhs = None, l
        rhs = rhs.strip()
        # print("cond : " + str(cond))
        if cond:
            cond = cond.strip()
        # print(rhs)
        d = match("_x | _y", rhs)
        if not d:
            lb = True
            d = match("_x | 0 <= _y", rhs)
        else:
            lb = False


        eq = d["_x"]
        comp = d['_y']
        var = comp
        # print( "... {}".format( (var, lb))  )

        conditions.append(cond)
        equations.append(eq)
        variables.append(var)
        complementarities.append(lb)

    return [conditions, equations, variables, complementarities]


def eval_ast(expr, d):
    import ast

    # a = ast.Expr(value=expr)
    a = ast.Expression(body=expr)
    a =     ast.fix_missing_locations(a)
    cc = compile(a, '<string>', 'eval')
    res = eval(cc, d)
    return res


def deterministic_simul(model, sol, ts_ind=0):
    import numpy
    variables = set([v.value.id for v in model.variables_ast])
    series = [ get_ts(model.etree, sol, v, ts_ind) for v in variables ]
    import pandas
    df = pandas.DataFrame( numpy.array(series).T, columns=variables)
    return df


class Model:
    __jacobian__ = None

    def __init__(self, etree, lines, calibration):

        [conditions, equations, variables, complementarities] = read_equations(lines)

        self.conditions_ast = conditions
        self.equations_ast = equations
        self.variables_ast = variables
        self.complementarities = complementarities
        self.calibration = calibration
        self.etree = etree
        pass


    def compute_jacobian(self, values={}):

        import time

        t1 = time.time()

        model = self
        etree = model.etree

        indexed_variables_str = set([va.value.id for va in model.variables_ast])
        indexed_variables = [LIVar(vs, etree, values) for vs in indexed_variables_str]

        context = {}
        for iv in indexed_variables:
            context[iv.name] = iv
        for c,v in self.calibration.items():
            context[c] = v

        t2 = time.time()
        import copy
        from .pattern import ReplaceExpectation

        def get_proba(etree,s):
            his = etree.history(s)
            transitions = [etree.probas[u,v] for u,v in zip(his[:-1],his[1:])]
            p = 1
            for pp in transitions:
                p *= pp
            return p

        context['s'] = etree.nodes[0]
        context['etree'] = etree
        context['S'] = lambda x: [(etree.probas[s, x], x) for x in etree.children(s)]
        context['P'] = lambda x: etree.parent(x)
        context['H'] = lambda x: etree.history(x)
        context['prob'] = lambda x: get_proba(etree, x)
        context['Ex'] = lambda fun, ss: sum([p*fun(x) for p,x in ss])
        context['Sum'] = lambda fun, ss: sum([fun(x) for x in ss])

        tsf_equations = [ReplaceExpectation().visit( copy.copy(ee) ) for ee in model.equations_ast]

        conditions_fun = [(eval("lambda s: " + l, context) if l is not None else None) for l in model.conditions_ast]


        # print(len(tsf_equations))
        t3 = time.time()

        full_equations = []
        full_variables = []
        full_complementarities = []

        # this will be quicker if I precompute substitutions.
        for n,eq in enumerate(tsf_equations):
            for s in etree.nodes:
                context['s'] = s
                context['t'] = len(s)-1
                if (conditions_fun[n] is None) or conditions_fun[n](s):
                    try:
                        var = eval_ast( model.variables_ast[n], context )
                    except Exception as e:
                        import ast
                        print( ast.dump(  model.variables_ast[n] ) )
                        # print()
                        raise e
                    full_variables.append(var.name)
                    eq_ast = tsf_equations[n]

                    eq = eval_ast(eq_ast, context)
                    # print( to_source( eq ) )
                    full_equations.append(  eq )
                    full_complementarities.append( model.complementarities[n] )
        t4 = time.time()

        P = len(full_equations)
        import numpy
        res = numpy.zeros(P)
        jac = numpy.zeros((P,P))
        lb = numpy.zeros(P)

        for p in range(P):
            eq = full_equations[p]
            res[p] = -eq.c[1]
            if full_complementarities[p]:
                lb[p] = 0
            else:
                lb[p] = -100000
            for k in eq.c.keys():
                if k != 1:
                    q = full_variables.index(k)
                    jac[p,q] = eq.c[k]

        t5 = time.time()

        # print('{} {} {} {}'.format(t2-t1,t3-t2,t4-t3,t5-t4))

        self.full_variables = full_variables
        self.full_equations = full_equations
        self.full_complementarities = full_complementarities

        return res, jac, lb

    def solve(self, verbose=False, jacs=None, linear=True, init_guess={}):

        if linear:

            model = self
            import numpy

            if jacs is None:
                constants,mat,lb = self.compute_jacobian()
            else:
                constants,mat,lb = jacs

            mm = numpy.array(mat).astype(dtype=float)
            v = numpy.array(constants).astype(dtype=float).flatten()
            lb = numpy.array(lb).astype(dtype=float).flatten()

            lb[lb<-1000] = -numpy.inf
            x0 = v*0 + 0.1
            ub = numpy.zeros_like(lb)+numpy.inf

            def fun(x):
                res = v-numpy.dot(mm, x)
                return res
            #
            def dfun(x):
                return -mm


            from dolo.numeric.extern.lmmcp import lmmcp
            res = lmmcp(fun, dfun, x0, lb, ub, verbose=verbose, options={'preprocess': True, 'eps2': 1e-10, 'presteps': 40})

            from collections import OrderedDict
            return OrderedDict(zip(map(str,model.full_variables), res))
        else:

            model = self
            import numpy

            constants, mat, lb = self.compute_jacobian(values=init_guess)   # warmup
            lb[lb < -1000] = -numpy.inf
            v = numpy.array(constants).astype(dtype=float).flatten()
            if len(init_guess)>0:
                x0 = numpy.array([init_guess[s] for s in self.full_variables])
            else:
                x0 = v*0 + 0.1
            ub = numpy.zeros_like(lb)+numpy.inf



            def resid(x):
                vals = {self.full_variables[i]: x[i] for i in range(len(x))}

                constants,mat,trash = self.compute_jacobian(values = vals)
                mm = numpy.array(mat).astype(dtype=float)
                v = numpy.array(constants).astype(dtype=float).flatten()
                # lb = numpy.array(lb).astype(dtype=float).flatten()


                # res = v-numpy.dot(mm, x)
                # return res
                return [v, -mm]

            def fun(x):
                return resid(x)[0]
            #
            def dfun(x):
                return resid(x)[1]


            from dolo.numeric.extern.lmmcp import lmmcp
            res = lmmcp(fun, dfun, x0, lb, ub, verbose=verbose, options={'preprocess': True, 'eps2': 1e-10})

            from collections import OrderedDict
            return OrderedDict(zip(map(str,model.full_variables), res))


if __name__ == '__main__':

    from .tree_import import import_tree_model
    lines = import_tree_model('models_tree.yaml','optimal')
    calibration = dict(
        beta = 0.8/(0.8+0.15),
        a = 0.8,
        c = 0.15,
        estar = 0.0,
        Rbar = 0.1,
        min_f = 0,
        kappa = 1.0,
        R = 0.00123
    )



    t0 = time.time()

    N = 20
        # etree = DeterministicTree(N)
    # for s in etree.nodes:
    #     etree.values[s] = 0.1
    etree = BranchTree(10, N, 0.5)
    for s in etree.nodes:
        if 1 in s:
            etree.values[s] = 0.1

    t1 = time.time()
    model = Model(etree, lines, calibration)
    t2 = time.time()
    res, jac, lb = model.compute_jacobian()
    t3 = time.time()
    res, jac, lb = model.compute_jacobian()
    t4 = time.time()
    sol = model.solve(verbose=False,jacs=[res,jac,lb])
    t5 = time.time()
    sim = deterministic_simul(model, sol)
    t6 = time.time()
    print("Total : {}".format(t5-t0))
    print('- construct the tree: {}'.format(t1-t0))
    print('- construct the model: {}'.format(t2-t1))
    print('- evaluate jacobien: {}'.format(t3-t2))
    print('- evaluate jacobien (2): {}'.format(t4-t3))
    print('- find solution: {}'.format(t5-t4))
    print('- simulate: {}'.format(t6-t5))

    t0 = time.time()

    t1 = time.time()
    model = Model(etree, lines, calibration)
    t2 = time.time()
    res, jac, lb = model.compute_jacobian()
    t3 = time.time()
    res, jac, lb = model.compute_jacobian()
    t4 = time.time()
    sol = model.solve(verbose=False,jacs=[res,jac,lb])
    t5 = time.time()
    sim = deterministic_simul(model, sol)
    t6 = time.time()
    print("Total : {}".format(t5-t0))
    print('- construct the tree: {}'.format(t1-t0))
    print('- construct the model: {}'.format(t2-t1))
    print('- evaluate jacobien: {}'.format(t3-t2))
    print('- evaluate jacobien (2): {}'.format(t4-t3))
    print('- find solution: {}'.format(t5-t4))
    print('- simulate: {}'.format(t6-t5))
