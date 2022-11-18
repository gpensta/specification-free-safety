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

    def draw(self, fig):
        plt.cla()
        plt.xlim(-1,1)
        plt.ylim(0, 1.2)
        plt.plot([0, self.x[0, 0]], [0, np.sqrt(1 - self.x[0, 0]**2)], color="black", linewidth=3, zorder=2)
        plt.scatter(self.x[0, 0], np.sqrt(1 - self.x[0, 0]**2), s=100, c='red', alpha=1., zorder=13)

    def move(self, u):
        self.x = np.array([[1., self.dt], [-self.dt * self.g/self.l, 1 -self.beta/self.m * self.dt]]) @ self.x + np.array([[0], [self.dt * u]])
    # def control(self):
    #     u = 10 * (-self.x[0, 0])
    #     return 7 * np.tanh(u)

    def control(self):
        a = -self.g/self.l
        b = -self.beta/self.m
        u = -(a+1) * self.x[0, 0] - (b+2) * self.x[1, 0]
        return 5 * np.tanh(u)
        


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
    u_list = []
    for i in range(20):
        if i < 5:
            u = 2 * np.random.rand() - 1 
        else:
            u = pendulum.control()
        u_list.append(u)
        pendulum.move(u)
        x1.append(pendulum.x[0, 0])
        x2.append(pendulum.x[1, 0])
        if (1 - pendulum.x[0, 0]**2) < 0:
            break
        pendulum.draw(fig)
        plt.savefig('results/pend.png')
        plt.pause(0.05)     
    #plots 
    plt.figure()
    plt.plot(x1)
    plt.savefig('results/x1.png')
    plt.figure()
    plt.plot(x2)
    plt.savefig('results/x2.png')
    plt.figure()
    plt.plot(u_list)
    plt.savefig('results/u.png')


def forward_reach(x, u):
    pass

def compute_weak_set():
    pass


if __name__ == '__main__':
    main(0.15)
    # draw_quiver()