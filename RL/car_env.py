import gym
from gym import spaces
from gym.utils import seeding

import numpy as np
import stable_baselines3
from stable_baselines3.common.env_checker import check_env
from stable_baselines3 import PPO
from stable_baselines3.common.env_util import make_vec_env
import matplotlib.pyplot as plt

class Car():
    def __init__(self, x1 = 0., x2 = 0.) -> None:
        self.x = np.array([[x1], [x2]])
        self.max_u = 1
        self.max_x = 10
        self.w = np.array([[0], [0]])
        self.dt = 1

    def draw(self, fig):
        plt.cla()
        plt.xlim(-1,1)
        plt.ylim(0, 1.2)
        plt.scatter(self.x[0, 0], self.x[1, 0], c='blue')
        plt.scatter(self.w[0, 0], self.w[1, 0], c='red')

    def move(self, u):
        self.x += self.dt + u.reshape((2,1))

class EnvCar(gym.Env):
    
    def __init__(self):
        self.agent = Car(0., 0.)
        self.step_nb = 0
        self.observation_space = spaces.Box(low=-self.agent.max_x, high=self.agent.max_x, shape=(2,), dtype=np.float32)
        self.action_space = spaces.Box(low=-self.agent.max_u, high=self.agent.max_u, shape=(2,), dtype=np.float32)
    
    def seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]
    
    def reset(self):
        self.step_nb = 0
        self.agent.x = 2 * self.agent.max_x * np.random.rand(2, 1) - self.agent.max_x
        return np.array([self.agent.x[0, 0], self.agent.x[1, 0]], dtype=np.float32)

    def supervisor(self):
        pass
    
    def step(self, action):
        self.step_nb += 1
        old_dist = np.linalg.norm(self.agent.x)
        self.agent.move(action)
        obs = np.array([self.agent.x[0, 0], self.agent.x[1, 0]], dtype=np.float32)
        dist = np.linalg.norm(self.agent.x)
        rew = old_dist - dist
        if dist < 0.1:
            return obs, 1., True, {}
        elif np.abs(self.agent.x[0, 0]) > self.agent.max_x or np.abs(self.agent.x[1, 0]) > self.agent.max_x:
            return obs, -1., True, {}
        elif self.step_nb > 25:
            return obs, -1., True, {}
        else:
            return obs, rew, False, {} 


def eval_model(model_path):
    model = PPO.load(model_path)
    env = EnvCar()
    success = 0
    for i in range(100):
        obs = env.reset()
        d = False
        while not d:
            action, _states = model.predict(obs)
            obs, reward, d, info = env.step(action)
        if reward > 0:
            success += 1
    return success / 100

if __name__ == '__main__':
    
    # env = EnvCar()
    # model = PPO("MlpPolicy", env, verbose=1)
    # model.learn(20000)
    # model.save("models/car")

    print(eval_model('/home/gwen/pendulum/RL/models/car.zip'))
    # model = PPO.load('models/inv_pend')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # fig = plt.figure()
    # x1 = (2 * np.random.rand() - 1.) * 0.3
    # pendulum = Pendulum(x1, 0)
    # pendulum.draw(fig)
    # x1 = []
    # x2 = []
    # u_list = []
    # for i in range(20):
    
    #     u, _ = model.predict(pendulum.x.flatten())
    #     u = u[0]
    #     u_list.append(u)
    #     pendulum.move(u)
    #     x1.append(pendulum.x[0, 0])
    #     x2.append(pendulum.x[1, 0])
    #     if (1 - pendulum.x[0, 0]**2) < 0:
    #         break
    #     pendulum.draw(fig)
    #     plt.savefig('results/pend.png')
    #     plt.pause(0.05)
    # plt.figure()
    # plt.plot(x1)
    # plt.savefig('results/x1.png')
    # plt.figure()
    # plt.plot(x2)
    # plt.savefig('results/x2.png')
    # plt.figure()
    # plt.plot(u_list)
    # plt.savefig('results/u.png')