# d = {}
# d["asd"] = 1
# print(d)
# print (range(1, 11))

import sys

print(sys.executable)
print(sys.version)


class MyClass:
    """docstring for MyClass"""

    def __init__(self, a):
        self.a = a

    def __eq__(self, other):
        print("MyClass.eq")
        return self.a == other.a

    def __hash__(self):
        print("MyClass.hash")
        return hash(self.a)

    def __repr__(self):
        return "{}".format(self.a)


print(list(range(10, 100 + 1, 10)))

o1 = MyClass(1)
o2 = MyClass(1)

print(repr(o1))

print(isinstance(o1, MyClass))
print(o1 == o2)

s = {MyClass(1), MyClass(2), MyClass(1)}
print(len(s))

objs = filter(lambda item: item.a == 1, s)

print(list(objs))

d = {}
for i in range(10):
    d[i] = MyClass(i)

print(range(0, 10, 1))

q = [1, 2, 3]
print (q)

print("==============NUMPY=====================")

import numpy as np
# print(np.sin(np.pi))

# print(np.cos(np.pi))

A = np.array([[1, 2, 3],
              [4, 5, 6],
              [7, 8, 9]])

# print(A[np.ix_([0, 2], [0, 1])])
print(A + 1)

print(np.insert(A, [0, 1, 3, 2], 5, axis=0))

print(list(range(10)))


v = {1, 2, 5, 3, 2}
print(v)
v = sorted(v)
print(v)

vec = np.zeros(3)
print(A.dot(vec).shape)

# noprimes = [j for i in range(2, 8) for j in range(i * 2, 50, i)]
# print(noprimes)
# primes = [x for x in range(2, 50) if x not in noprimes]
# print(primes)


# from scipy.integrate import quad as q
# # from scipy.integrate import fixed_quad as gq


# # def myfunc2(x, y):
# #     print("myfunc2: y={}".format(y))
# #     return x + y


# # def integr2(y):
# #     print("integr2: y={}".format(y))
# #     return q(myfunc2, 0, 1, y)
# #     # return q(myfunc2, 0, 1, args=[y])[0]


# # # q(myfunc2, 0, 1)
# # print(q(integr2, 0, 1))

# def f(x, y):
#     return x + y


# print(q2d(f, 0, 1, 0, 1))

# # print(gq(lambda x: (f(x[0], x[1]), 2), 0, 1))
