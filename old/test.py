import numpy as np

def simulate_window_observation(sequence_size, window_size, num_position, num_simulations):
    count_close = []
    mean_dist = []
    positions = list(range(1, sequence_size+1))
    for sim in range(num_simulations):
        # Sample random positions
        rand_pos = list(np.random.choice(positions, replace=False, size=num_position))
        rand_pos.sort()

        # Check if at least the difference is less than or equal to the window size
        close = 0
        all_dist = []
        pos0 = rand_pos[0]
        for i in rand_pos[1:]:
            dist = i - pos0
            all_dist.append(dist)

            if dist <= window_size:
                close += 1
            pos0 = i

        count_close.append(close)
        mean_dist.append(np.mean(all_dist))

    # Calculate the proportion of TRUE
    # proportion_true = count_true / num_simulations
    # proportion_false = 1 - proportion_true
    # return proportion_true, proportion_false
    return count_close, mean_dist

# Example with a 100-bp sequence and a 10-bp window
sequence_size = 250
window_size = 10
num_simulations = 10000
num_position = 10

count_close, mean_dist = simulate_window_observation(sequence_size, window_size, num_position, num_simulations)

print(f"Mean close: {np.mean(count_close):.4f}")
print(f"Mean dist: {np.mean(mean_dist):.4f}")


#print(f"Proportion of TRUE: {proportion_true:.4f}")
#print(f"Proportion of FALSE: {proportion_false:.4f}")

