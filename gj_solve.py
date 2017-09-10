# -*- coding: utf-8 -*-

def shape(M):
    return len(M),len(M[0]) if len(M) > 0 else 0

def matxRound(M, decPts=4):
    for row in M:
        for i,col in enumerate(row):
            row[i] = round(col, decPts)

def transpose(M):
    if len(M) == 0 or len(M[0]) == 0:
        return M
    row_n = len(M)
    col_n = len(M[0])
    MT = [[0] * row_n for i in xrange(col_n)]
    for i in xrange(row_n):
        for j in xrange(col_n):
            MT[j][i] = M[i][j]
    return MT

def matxMultiply(A, B):
    row_a = len(A)
    if row_a == 0:
        return None
    col_a = len(A[0])
    if col_a == 0:
        return None
    
    row_b = len(B)
    if row_b == 0:
        return None
    col_b = len(B[0])
    if col_b == 0:
        return None
    
    if col_a != row_b:
        return None
    
    M = [[0] * col_b for i in xrange(row_a)]
    
    for i in xrange(row_a):
        for j in xrange(col_b):
            M[i][j] = sum([A[i][k] * B[k][j] for k in xrange(col_a)])
    
    return M

def augmentMatrix(A, b):
    if len(A) == 0:
        return A
    AM = [[0] * (1+len(A[0])) for i in xrange(len(A))]
    for i in xrange(len(A)):
        for j in xrange(len(A[0])):
            AM[i][j] = A[i][j]
    for i in xrange(len(A)):
        AM[i][len(A[0])] = b[i][0]

    return AM

# TODO r1 <---> r2
# 直接修改参数矩阵，无返回值
def swapRows(M, r1, r2):
    M[r1], M[r2] = M[r2], M[r1]

# TODO r1 <--- r1 * scale， scale!=0
# 直接修改参数矩阵，无返回值
def scaleRow(M, r, scale):
    if scale == 0:
        raise ValueError
    for i,v in enumerate(M[r]):
        M[r][i] = v * scale

# TODO r1 <--- r1 + r2*scale
# 直接修改参数矩阵，无返回值
def addScaledRow(M, r1, r2, scale):
    if scale == 0:
        return
    for i,v in enumerate(M[r1]):
        M[r1][i] = v + M[r2][i]*scale

""" Gaussian Jordan 方法求解 Ax = b.
    参数
        A: 方阵 
        b: 列向量
        decPts: 四舍五入位数，默认为4
        epsilon: 判读是否为0的阈值，默认 1.0e-16
        
    返回列向量 x 使得 Ax = b 
    返回None，如果 A，b 高度不同
    返回None，如果 A 为奇异矩阵
"""
#     对于Ab的每一列（最后一列除外）
#         当前列为列c
#         寻找列c中 对角线以及对角线以下所有元素（行 c~N）的绝对值的最大值
#         如果绝对值最大值为0
#             那么A为奇异矩阵，返回None （请在问题2.4中证明该命题）
#         否则
#             使用第一个行变换，将绝对值最大值所在行交换到对角线元素所在行（行c） 
#             使用第二个行变换，将列c的对角线元素缩放为1
#             多次使用第三个行变换，将列c的其他元素消为0
#             
# 步骤4 返回Ab的最后一列
def gj_Solve(A, b, decPts=4, epsilon = 1.0e-16):
    if len(A) != len(b):
        return None

    if len(A) == 0 or len(A[0]) == 0:
        return None

    import pprint
    pp = pprint.PrettyPrinter(indent=1,width=20)

    rows = len(A)
    cols = len(A[0])

    # 非方阵
    if rows != cols:
        return None

    AM = augmentMatrix(A, b)
    for c in xrange(cols):
        # find max abs
        max_abs_row = c
        max_abs = abs(AM[c][c])
        for i in xrange(c+1, rows):
            if abs(AM[i][c]) > max_abs:
                max_abs = abs(AM[i][c])
                max_abs_row = i
        # max abs 为0，是奇异矩阵，返回None
        if (max_abs - 0) <= epsilon:
            return None
        swapRows(AM, max_abs_row, c)

        # scale to 1
        scaleRow(AM, c, 1./AM[c][c])

        # clear other rows
        for i in xrange(rows):
            if i == c:
                continue
            addScaledRow(AM, i, c, -AM[i][c])   
    
    matxRound(AM, decPts)

    return [[AM[i][cols]] for i in xrange(rows)]


def matxLossMean(A, B):
    if len(A) != len(B):
        return None
    if len(A[0]) != len(B[0]):
        return None

    loss = 0.
    for i in xrange(len(A)):
        for j in xrange(len(A[0])):
            loss += (A[i][j] - B[i][j]) ** 2
    
    return loss / (len(A)*len(A[0]))

def test_case():
    import pprint
    pp = pprint.PrettyPrinter(indent=1,width=20)

    # TODO 构造 矩阵A，列向量b，其中 A 为奇异矩阵
    A = [[1,2,3],
         [0,0,5],
         [0,0,1]]
    b = [[4],
         [5],
         [6]]
    x = gj_Solve(A, b)
    print 'singular matrix passed?', x == None

    # TODO 构造 矩阵A，列向量b，其中 A 为非奇异矩阵
    A = [[5.262, 2.739, -9.878],
         [5.111, 6.358, 7.638],
         [2.016, -9.924, -1.367]]
    b = [[-3.441],
         [-2.152],
         [-9.278]]

    # TODO 求解 x 使得 Ax = b
    x = gj_Solve(A, b)
    print 'x='
    pp.pprint(x)

    # TODO 计算 Ax
    Ax = matxMultiply(A, x)
    print 'Ax='
    pp.pprint(Ax)

    # TODO 比较 Ax 与 b
    print 'Ax == b?', matxLossMean(Ax, b) <= 1.0e-6



'''
参数：(x,y) 二元组列表
返回：m，b
'''
def linearRegression(points):
    # XTXh=XTY
    X = []
    Y = []
    for x,y in points:
        X.append([x, 1.])
        Y.append([y])
    h = gj_Solve(matxMultiply(transpose(X), X), matxMultiply(transpose(X), Y))
    if h is None:
        return 0,0
    else:
        return h[0][0], h[1][0]

# TODO 构造线性函数
def f(x):
    """$y = 4x + 3$"""
    return (4*x + 3)

def test_case2():
    # TODO 构造 100 个线性函数上的点，加上适当的高斯噪音
    import random
    points = []
    gauss_sigma = 1
    for i in xrange(100):
        x = random.randint(0, 100)
        points.append((x+random.gauss(0, gauss_sigma), f(x)+random.gauss(0, gauss_sigma)))

    #TODO 对这100个点进行线性回归，将线性回归得到的函数和原线性函数比较
    h = linearRegression(points)
    print "linearRegression(coef, intercept):", h

    import matplotlib.pyplot as plt
    import numpy as np
    px = np.arange(0.0, 100.0, 0.01)
    py = f(px)
    plt.title(f.__doc__, fontsize=20)
    plt.scatter([i[0] for i in points], [i[1] for i in points])
    plt.plot(px, py, color='r')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.text(10, 5, 'y={}x + {}'.format(h[0], h[1]))
    plt.show()

test_case()
test_case2()

