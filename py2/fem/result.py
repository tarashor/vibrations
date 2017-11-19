import numpy as np
from . import finiteelements as fe

class Result:
    def __init__(self, freq, u1, u2, u3, mesh):
        self.freq = freq
        self.u1 = u1
        self.u2 = u2
        self.u3 = u3
        self.mesh = mesh

    def get_displacement(self, x1, x2, x3, time):
        element = self.mesh.get_element(x1, x2)

        u_res = np.zeros((3))

        h_e = lin_aprox_matrix(node, x1, x2, x3)
        u_res += u_node * f


        return u_res * self.fi(time)

    def get_u_in_node(self, node):
        u = np.zeros((3))
        u[0] = self.u1[node.index]
        u[1] = self.u2[node.index]
        u[2] = self.u3[node.index]
        return u

    def fi(self, time):
        return np.sin(self.freq * time)


class Result:
    def __init__(self, lam, vec, mesh, model):
        self.lam = lam
        self.vec = vec
        self.mesh = mesh
        self.model = model

    def get_results_count(self):
        return len(self.lam)

    def get_nodes(self):
        return self.mesh.nodes

    def get_result(self, i):
        return np.sqrt(self.lam[i]), self.vec[:, i][0:self.mesh.nodes_count()], self.vec[:, i][self.mesh.nodes_count():2 * self.mesh.nodes_count()], self.mesh.nodes

    def get_result_min(self):
        return self.get_result(0)


def get_strain(self, freq_index, alpha1, alpha2):
        res = self.vec[:, freq_index]
        normalize(res)
        element = self.mesh.get_element(alpha1, alpha2)
        geometry = self.model.geometry

        B = deriv_to_grad(alpha1, alpha2, geometry)
        I_e = ksiteta_to_alpha_matrix(element)

        ksi = 2 * alpha1 / element.width() - (element.top_left.x + element.top_right.x) / element.width()
        teta = 2 * alpha2 / element.height() - (element.top_left.y + element.bottom_left.y) / element.height()

        H = lin_aprox_matrix(ksi, teta)

        u_e = np.zeros(8)

        for i in range(8):
            i_g = map_local_to_global_matrix_index(i, element, res.shape[0])
            u_e[i] = res[i_g]

        grad_u = B.dot(I_e).dot(H).dot(u_e)

        E = grad_to_strain_linear_matrix()
        E_NL = grad_to_strain_nonlinear_matrix(alpha1, alpha2, geometry, grad_u)
        return (E+0.5*E_NL).dot(grad_u)


print(np.ones((2,3)).dot(np.zeros((3))))
