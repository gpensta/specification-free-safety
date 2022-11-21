import gym
from gym import spaces
from gym.utils import seeding
from discrete_pendulum import Pendulum
import numpy as np
import stable_baselines3
from stable_baselines3.common.env_checker import check_env
from stable_baselines3 import PPO
from stable_baselines3.common.env_util import make_vec_env


class EnvPendulum(gym.Env):
    
    def __init__(self):
        self.agent = Pendulum(0., 0.)
        self.step_nb = 0
        self.observation_space = spaces.Box(low=-1, high=1, shape=(2,), dtype=np.float32)
        self.action_space = spaces.Box(low=-5., high=5., shape=(1,), dtype=np.float32)
    
    def seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]
    
    def reset(self):
        x1 =  (2 * np.random.rand() - 1.) * 0.3
        self.agent.x[0, 0] = x1
        self.agent.x[1, 0] = 0.
        return np.array([self.agent.x[0, 0], self.agent.x[1, 0]], dtype=np.float32)
    
    def step(self, action):
        self.step_nb += 1
        old_dist = np.abs(self.agent.x[0, 0])
        self.agent.move(action[0])
        obs = np.array([self.agent.x[0, 0], self.agent.x[1, 0]], dtype=np.float32)
        dist = np.abs(self.agent.x[0, 0])
        rew = old_dist - dist
        if np.abs(self.agent.x[0, 0]) < 0.1 and np.abs(self.agent.x[1, 0]) < 0.1:
            return obs, 1., True, {}
        elif np.abs(self.agent.x[0, 0]) > 0.45:
            return obs, -1., True, {}
        else:
            return obs, rew, False, {} 


def eval_model(model_path):
    model = PPO.load(model_path)
    env = EnvPendulum()
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
    env = EnvPendulum()
    model = PPO("MlpPolicy", env, verbose=1)
    model.learn(100000)
    model.save("models/inv_pend")