import numpy as np
import matplotlib.pyplot as plt 

class Pendulum():
    def __init__(self, x1, x2) -> None:
        self.x = np.array([[x1], [x2]])
        self.l = 1
        self.g = -10
        self.beta = 1
        self.m = 1
        self.dt = 0.1

    def f(self, x, u):
        A = np.array([[0, 1], [-self.g/self.l, -self.beta/self.m]])
        return A @ x + u * np.array([[0], [1]]) 

    def draw(self, fig):
        plt.cla()
        plt.xlim(-1,1)
        plt.ylim(0, 1.2)
        plt.plot([0, self.x[0, 0]], [0, np.sqrt(1 - self.x[0, 0]**2)])

    def move(self, u):
        self.x += self.dt * self.f(self.x, u)

    # def control(self):
    #     u = 10 * (-self.x[0, 0])
    #     return 7 * np.tanh(u)

    def control(self):
        a = -self.g/self.l
        b = -self.beta/self.m
        u = -(a+1) * self.x[0, 0] - (b+2) * self.x[1, 0]
        return 2 * np.tanh(u)


def draw_quiver():
    x1_list = []
    x2_list = []
    x1_direction_list = []
    x2_direction_list = []
    colors = []
    for x1 in np.linspace(-1. , 1., 30):
        for x2 in np.linspace(-5. , 5, 30):
            x1_list.append(x1)
            x2_list.append(x2)
            pend = Pendulum(x1, x2)
            u = pend.control()
            pend.move(u)
            norm = np.linalg.norm(np.array([(pend.x[0, 0] - x1), (pend.x[1, 0] - x2)]))
            colors.append(norm)
            x1_direction_list.append((pend.x[0, 0] - x1) / norm)# / np.linalg.norm((pend.x[0, 0] - x1)))
            x2_direction_list.append((pend.x[1, 0] - x2)/ norm)# / np.linalg.norm(pend.x[1, 0] - x2))
    fig, ax = plt.subplots(figsize = (12, 7))
    ax.quiver(x1_list, x2_list, x1_direction_list, x2_direction_list,
         colors, scale = 30)
    plt.savefig('quiver.png')
 
def main(x1):
    fig = plt.figure()
    pendulum = Pendulum(x1, 0)
    pendulum.draw(fig)
    x1 = []
    x2 = []
    for i in range(100):
        u = pendulum.control()
        pendulum.move(u)
        x1.append(pendulum.x[0, 0])
        x2.append(pendulum.x[1, 0])
        if (1 - pendulum.x[0, 0]**2) < 0:
            break
        pendulum.draw(fig)
        plt.savefig('pend.png')
        print(u)
        plt.pause(0.05)     
    plt.figure()
    plt.plot(x1)
    plt.savefig('x1.png')
    plt.figure()
    plt.plot(x2)
    plt.savefig('x2.png')


def forward_reach(x, u):
    pass

def compute_weak_set():
    pass


if __name__ == '__main__':
    # main(0.15)
    draw_quiver()