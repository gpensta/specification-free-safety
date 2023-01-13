import numpy as np

# Generate some toy data to fit the GP to
np.random.seed(42)
X = np.random.uniform(-3, 3, (20, 1))
Y = np.sin(X) + np.random.randn(20, 1) * 0.1

# Define the mean and covariance functions
def mean_function(x):
    return np.zeros(x.shape[0])

def covariance_function(x, x_prime):
    return np.exp(-0.5 * np.sum((x - x_prime) ** 2, axis=1))

# Make predictions at new points
x_new = np.linspace(-3, 3, 20).reshape(-1, 1)

# Compute the mean and covariance of the posterior distribution
K = covariance_function(X, X) + 1e-8 * np.eye(X.shape[0])
K_s = covariance_function(X, x_new)
K_ss = covariance_function(x_new, x_new)

mean_post = mean_function(x_new) + K_s.T.dot(np.linalg.inv(K)).dot(Y - mean_function(X))
cov_post = K_ss - K_s.T.dot(np.linalg.inv(K)).dot(K_s)

# Draw samples from the posterior distribution
samples = np.random.multivariate_normal(mean_post.flatten(), cov_post, 3)

# Plot the results
import matplotlib.pyplot as plt
plt.figure()
plt.plot(x_new, mean_post, 'k', lw=2, zorder=9)
plt.fill_between(x_new[:, 0], mean_post[:, 0] - 1.96 * np.sqrt(cov_post.diagonal()), mean_post[:, 0] + 1.96 * np.sqrt(cov_post.diagonal()), alpha=0.2, color='k')
plt.scatter(X, Y, c='r', s=50, zorder=10, edgecolors=(0, 0, 0))
for i in range(samples.shape[0]):
    plt.plot(x_new, samples[i, :], color='r', linewidth=1)
plt.title("GP regression")
plt.savefig("gp.png")
