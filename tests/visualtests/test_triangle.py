import numpy as np
import matplotlib.pyplot as plt

class Triangle:

    def __init__(self, ver, idx, orient, plot_triangle=False):
        orient = np.array(orient)
        A = ver[idx[0]]
        B = ver[idx[1]]
        C = ver[idx[2]]

        v1 = B - A
        v2 = C - A
        normal = np.cross(v1, v2)
        normal = normal / np.linalg.norm(normal)

        orient = orient - np.dot(orient, normal) * normal
        orient = orient / np.linalg.norm(orient)

        orthogonal = np.cross(normal, orient)
        s0 = 0
        s1 = np.dot(v1, orthogonal)
        s2 = np.dot(v2, orthogonal)
        #print('%.4f %.4f' % (s1, s2))

        i0 = np.argmin([s0, s1, s2])
        i1 = np.argmax([s0, s1, s2])
        i2 = 3-i0-i1
        if (i1-i0) % 3 == 2:
            i0, i1 = i1, i0
        #print(i0, i1, i2)

        idx_origin, idx_base, idx_apex = idx[i0], idx[i1], idx[i2]
        x_origin = ver[idx_origin]
        x_base = ver[idx_base]
        x_apex = ver[idx_apex]
        
        OB = x_base - x_origin
        OA = x_apex - x_origin

        M = [[np.dot(orient, OB), 1],
             [np.dot(OB, OB), np.dot(orient, OB)]]
        s_crit, height = np.linalg.solve(M, [np.dot(orient, OA), np.dot(OA, OB)])
        assert 0 < s_crit < 1
        if height < 0:
            orient = -orient
            height = -height
        #print('%.4f %.4f' % (s_crit, height))
        #print('%.4f %.4f %.4f' % (orient[0], orient[1], orient[2]))
        if plot_triangle:
            plt.plot(x_origin[0], x_origin[1], 'ro')
            plt.plot(x_base[0], x_base[1], 'bo')
            x = np.row_stack((x_origin, x_origin+OB*s_crit, x_origin+OB*s_crit+height*orient))
            plt.plot(x[:, 0], x[:, 1], c='g')
            assert np.allclose(x[2], x_apex)
            assert np.dot(normal, np.cross(OB, OA)) > 0
            x = np.row_stack((A, B, C, A))
            plt.plot(x[:, 0], x[:, 1], ls='--', c=(0.8, 0.8, 0.8))
            G = (A+B+C)/3
            plt.plot([G[0]-orient[0], G[0]+orient[0]], [G[1]-orient[1], G[1]+orient[1]])
            #plt.plot([G[0]-orthogonal[0], G[0]+orthogonal[0]], [G[1]-orthogonal[1], G[1]+orthogonal[1]])
            plt.text(A[0], A[1], 'A')
            plt.text(B[0], B[1], 'B')
            plt.text(C[0], C[1], 'C')
            plt.axis('equal')
            plt.show()

        self.ver = ver
        self.idx_origin = idx_origin
        self.idx_base = idx_base
        self.idx_apex = idx_apex
        self.origin = x_origin
        self.base_vector = OB
        self.orient = orient
        self.s_crit = s_crit
        self.height = height
    
    def segment_length(self, s):
        if s <= self.s_crit:
            return self.height*s/self.s_crit
        else:
            return self.height*(1-s)/(1-self.s_crit)
    
    def segment(self, s):
        P = self.origin + s*self.base_vector
        Q = P + self.segment_length(s) * self.orient
        return P, Q
    
    def plot_tri(self):
        x = np.row_stack((
            self.ver[self.idx_origin],
            self.ver[self.idx_base],
            self.ver[self.idx_apex],
            self.ver[self.idx_origin],
        ))
        plt.plot(x[:, 0], x[:, 1], ls='--', c=(0.8, 0.8, 0.8))
        plt.text(x[0, 0], x[0, 1], 'O')
        plt.text(x[1, 0], x[1, 1], 'B')
        plt.text(x[2, 0], x[2, 1], 'A')
    
    def plot_segments(self, n=10):
        for s in np.linspace(0, 1, n):
            P, Q = self.segment(s)
            assert np.abs(np.linalg.norm(P-Q)- self.segment_length(s)) < 1e-6
            plt.plot([P[0], Q[0]], [P[1], Q[1]], 'r')

    
def generate_data():
    ver = np.random.random(size=(100, 3))
    ver.tofile('data_ver.txt', sep=' ')
    itri = np.random.randint(ver.shape[0], size=(100, 3))
    degenerate = (itri[:, 0] == itri[:, 1]) | (itri[:, 0] == itri[:, 2]) | (itri[:, 2] == itri[:, 1])
    itri = itri[~degenerate]
    itri = itri[:50]
    itri.tofile('data_itri.txt', sep=' ')
    np.random.random(size=(50, 3)).tofile('data_orient.txt', sep=' ')


def test_triangle():
    ver = np.fromfile('data_ver.txt', sep=' ').reshape((-1, 3))
    tri = np.fromfile('data_itri.txt', sep=' ').reshape((-1, 3)).astype(int)
    orient = np.fromfile('data_orient.txt', sep=' ').reshape((-1, 3))
    nt = tri.shape[0]
    for i in range(nt):
        T = Triangle(ver, tri[i], orient[i])
        #T.plot_tri()
        #T.plot_segments(50)
        #plt.show()
        P, Q = T.segment(0.4)
        print('%.4f %.4f %.4f' % (Q[0], Q[1], Q[2]))


def test_mask_sphere():
    x = np.array([[2, 1, 0],
                  [5, 1, 0],
                  [3, 3, 0]])
    x = np.vstack((x, x[:1]))
    c = np.array([4, 0, 0])
    R = 1.25
    plt.plot(x[:, 0], x[:, 1])
    t = np.linspace(0, 2*np.pi, 100)
    plt.plot(c[0]+R*np.cos(t), c[1]+R*np.sin(t))
    plt.axis('equal')
    plt.show()


def test_mask_cylinder():
    x = np.array([[2, 1, 0],
                  [5, 1, 0],
                  [3, 3, 0]])
    x = np.vstack((x, x[:1]))
    #PQ = np.array([[2, 0, 0], [5, 3, 0]])
    PQ = np.array([[2, 0, 0], [3.5, 1.5, 0]])
    #PQ = np.array([[2, 2, 0], [4.5, 2, 0]])
    #PQ = np.array([[2.5, 1.5, 0], [3.5, 2, 0]])
    R = 0.5
    plt.plot(x[:, 0], x[:, 1])

    t = np.linspace(0, 2*np.pi, 100)
    plt.plot(PQ[:, 0], PQ[:, 1])
    plt.plot(PQ[0, 0]+R*np.cos(t), PQ[0, 1]+R*np.sin(t))
    plt.plot(PQ[1, 0]+R*np.cos(t), PQ[1, 1]+R*np.sin(t))
    d = PQ[1] - PQ[0]
    d = d / np.linalg.norm(d)
    n = np.array([d[1], -d[0], 0])
    plt.plot(PQ[:, 0]+R*n[0], PQ[:, 1]+R*n[1])
    plt.plot(PQ[:, 0]-R*n[0], PQ[:, 1]-R*n[1])

    plt.axis('equal')
    plt.show()


def test_linear_constraint():
    x = np.array([[2, 1, 0],
                  [5, 1, 0],
                  [3, 3, 0]])
    x = np.vstack((x, x[:1]))
    PQ = np.array([[3, 0.5, 0], [4, 1.5, 0]])
    plt.plot(x[:, 0], x[:, 1])
    plt.plot(PQ[:, 0], PQ[:, 1], 'o-')

    a, b, c = 1.50000, 1.00000, -0.25000
    s = np.linspace(0, 1, 10)
    t = -(a*s+c)/b
    plt.plot(2+3*s, 1+2*t)
    t = -(a*s+c-1)/b
    plt.plot(2+3*s, 1+2*t)
    plt.axis('equal')
    plt.show()


#test_linear_constraint()
test_mask_cylinder()