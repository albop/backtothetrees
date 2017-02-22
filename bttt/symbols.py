class Symbol:

    def __init__(self, name=None, subscripts=None, date=None):
        if date is None:
            date = 0
        if subscripts is None:
            subscripts = []
        assert(isinstance(date,int))
        for s in subscripts:
            assert(isinstance(s,str))
            assert('__' not in s)
        assert(name.isidentifier())
        self.subscripts = tuple( subscripts )
        self.date = date
        self.name = name

    def __str__(self):
        # sep = 'ⵂ'
        # close_lhs = 'ⵂ'
        # close_rhs = 'ⵂ'
        sep = '__'
        close_lhs = '__'
        close_rhs = ''
        # close_lhs = '᚜'
        # close_rhs = '᚛'
        if self.date:
            argss = str.join(sep, [write_time_subscript(self.date)] + list(self.subscripts))
        else:
            argss =  str.join(sep, list(self.subscripts))
        s = self.name + close_lhs + argss + close_rhs
        return s

def write_time_subscript(s):
    # neg = '֊'
    neg = 'm'
    if isinstance(s, int):
        if s>=0:
            return str(s)
        elif s<0:
            return neg + str(-s)
    else:
        raise Exception("Don't know what to do with that")

class LinearExpression:

    def __init__(self, c=None, name=None, values={}):
        if name is not None:
            c = {name: 1.0}
        if c is None:
            c = dict()
        if 1 not in c:
            val = values.get(name, 0.0)
            c[1] = val
        self.c = c
        self.name = name

    def __add__(self, d):
        kk = self.c.copy() # ?
        if isinstance(d, LinearExpression):
            for (k,v) in d.c.items():
                if k in kk:
                    kk[k] = kk[k] + v
                else:
                    kk[k] = v
        else:
            kk[1] = kk[1] + d
        return LinearExpression(c=kk)

    __radd__ = __add__


    def __sub__(self, d):
        kk = self.c.copy() # ?
        if isinstance(d, LinearExpression):
            for (k,v) in d.c.items():
                if k in kk:
                    kk[k] = kk[k] - v
                else:
                    kk[k] = -v
        else:
            kk[1] = kk[1] - d
        return LinearExpression(c=kk)

    def __rsub__(self, d):
        return self.__neg__() + d

    def __neg__(self):
        kk = self.c.copy() # ?
        for k in kk.keys():
            kk[k] = -kk[k]
        return LinearExpression(c=kk)

    def __mul__(self, v):
        if isinstance(v, LinearExpression):
            a = self
            b = v
            c_a = a.c[1]
            c_b = b.c[1]
            res = c_a*b + c_b*a
            res.c[1] = c_a*c_b
            return res
        else:
            kk = self.c.copy() # ?
            for k in kk.keys():
                kk[k] *= v
            return LinearExpression(c=kk)

    __rmul__ = __mul__

    def __truediv__(self, v):
        kk = self.c.copy() # ?
        for k in kk.keys():
            kk[k] /= v

        return LinearExpression(c=kk)

    def __repr__(self):
        return str(self.c)



class LIVar:
    # indexed symbol which returns a linear expression

    def __init__(self, name, etree, values):
        self.name = name
        self.etree = etree
        self.values = values

    def __getitem__(self, x):
        ind = self.etree.nodes.index(x)   # why use the index ?
        name = '{}_{}'.format(self.name,ind)
        return LinearExpression(name=name, values=self.values)
