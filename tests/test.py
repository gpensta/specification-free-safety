import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def clip(x, inf, sup):
    if x < inf:
        x = inf
    elif x > sup:
        x = sup 
    return x

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

    def pd_control(self, sup_u):
        a = -self.g/self.l
        b = -self.beta/self.m
        u = -(a+1) * self.x[0, 0] - (b+2) * self.x[1, 0]
        return clip(u, -sup_u, sup_u)

def draw_streamlines():
    n = 500
    pend = Pendulum(0, 0)
    X, Y = np.meshgrid(np.linspace(-1, 1, n), np.linspace(-5, 5, n))
    U = np.ones((n, n))
    V = np.ones((n, n))
    for i in range(n):
        for j in range(n):
            x1 = X[i, j]
            x2 = Y[i, j]
            pend.x =np.array([[x1], [x2]])
            u = pend.pd_control(5.)
            pend.move(u)
            U[i, j] = (pend.x[0, 0] - x1)
            V[i, j] = (pend.x[1, 0] - x2)
    fig, ax = plt.subplots(figsize = (12, 7))
    # ax.quiver(x1_list, x2_list, x1_direction_list, x2_direction_list,
    ax.streamplot(X, Y, U, V, density=1.)
    plt.savefig('tests/streamline.png')

def draw_quiver():
    x1_list = []
    x2_list = []
    x1_direction_list = []
    x2_direction_list = []
    colors = []
    ntraj = 10
    traj = np.ones((2, ntraj))
    pend = Pendulum(0, 0)
    for x1 in np.linspace(-1. , 1., 30):
        for x2 in np.linspace(-5. , 5, 30):
            x1_list.append(x1)
            x2_list.append(x2)
            pend.x =np.array([[x1], [x2]])
            u = pend.pd_control(10)
            pend.move(u)
            norm = np.linalg.norm(np.array([(pend.x[0, 0] - x1), (pend.x[1, 0] - x2)]))
            colors.append(norm)
            x1_direction_list.append((pend.x[0, 0] - x1) / norm)
            x2_direction_list.append((pend.x[1, 0] - x2)/ norm)
    pend.x =np.array([[-0.25], [-2.]])
    for i in range(ntraj):
        traj[0, i] = pend.x[0, 0]
        traj[1, i] = pend.x[1, 0]
        u = pend.pd_control(10)
        pend.move(u)

    fig, ax = plt.subplots(figsize = (12, 7))
    ax.quiver(x1_list, x2_list, x1_direction_list, x2_direction_list,
         colors, scale = 30)
    ax.plot(traj[0, :-1], traj[1, :-1])
    
    plt.savefig('tests/state_space.png')

# draw_quiver()
draw_streamlines()

# time_horizon = 60
# max_u = 5


# n = 10
# heat_map = np.ones((n, n))

# X1 = np.linspace(-1, 1, n)
# X2 = np.linspace(-5, 5, n)

# pend = Pendulum(0, 0)

# for i in range(n):
#     for j in range(n):
#         x2 = X2[-(i+1)]
#         x1 = X1[j]
#         pend.x = np.array([[x1], [x2]])
#         for t in range(time_horizon):
#             u = pend.pd_control(5)
#             pend.move(u)
#         dist = np.abs(np.linalg.norm(pend.x))
#         if dist < np.linalg.norm(np.array([[x1], [x2]])):
#             dist = 0
#         else:
#             dist = 1
#         heat_map[i, j] = dist

# # print(heat_map)
# plt.figure()
# # plt.imshow(heat_map, cmap='hot', interpolation='nearest')
# # from matplotlib.colors import LinearSegmentedColormap

# # cmap_reds = plt.get_cmap('Reds')
# # num_colors = 15
# # colors = ['blue'] + [cmap_reds(i / num_colors) for i in range(2, num_colors)]
# # cmap = LinearSegmentedColormap.from_list('', colors, num_colors)
# # ax = sns.heatmap(heat_map, cmap=cmap, vmin=0, vmax=np.max(heat_map), square=True, cbar=False, annot=True)
# # cbar = plt.colorbar(ax.collections[0], ticks=range(num_colors + 1))

# plt.title
# ax = sns.heatmap(heat_map)
# ax.set(xticklabels=[])
# ax.set(yticklabels=[])
# ax.tick_params(bottom=False, left=False)
# plt.xlabel('x1 = [-1, 1]')
# plt.ylabel('x2 = [-5, 5]')



# # plt.contourf(heat_map, levels=[2,3])
# # sns.kdeplot(heat_map)
# plt.savefig("tests/heat_map.png")




