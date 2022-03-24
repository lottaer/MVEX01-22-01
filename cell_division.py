import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn; sn.set_theme()
import time
import pprint

class Ticker:
    def __init__(self):  # methods. The ticker is giving the ID
        self.counter = 0

    def get(self):
        self.counter += 1
        return self.counter

    def reset(self):
        self.counter = 0

def divide_cell(population, ticker, cell, time_units, p, k, n):
    bio_age = [cell[0]]

    cell[3] = ticker.get()  # [0, 0, 0, 0] cell type(bio), age, mother type, ID

    for i in range(time_units):

        # accumulate n new proteins
        proteins = cell[0]
        proteins += n

        # calculate how many proteins to transfer
        n_of_proteins_to_transfer = np.random.binomial(proteins, p)

        # divide (create new) cell
        new_cell = [n_of_proteins_to_transfer, 0, proteins-n, 0]  # cell type(bio), age, mother type, ID
        proteins -= n_of_proteins_to_transfer  # in mother cell
        cell[1] += 1  # increase chronological age
        if new_cell[0] <= k:
            # divide_cell(population, ticker, new_cell, time_units-cell[1], p, k, n)
            pass
        if proteins <= k:
            cell[0] = proteins
            bio_age.append(proteins)
        else:
            cell[0] = proteins
            # print(cell, 'Dead at age', cell[1])
            break
    population.append(cell)
    result = {'bio_age': bio_age,
              'age': len(bio_age)-1}
    return result


start = time.process_time()  # measuring the speed of the code

trials = 1000
time_units = 100
mothers = np.zeros((trials, time_units+1))  # creating zero matrix
mothers_lifespan = np.zeros(trials)  # for making average lifespan of the mothers
ps = [0.2, 0.5, 0.8]  # 3 different iterations
k = 3
n = 2

mean_mothers = np.zeros((len(ps), time_units+1))
ages_mothers = np.zeros((len(ps), time_units+1))

for iteration, p in enumerate(ps):
    print('iteration:', iteration)
    print('p:', p)
    for trial in range(0,trials):
        ticker = Ticker()
        population = []
        result = divide_cell(population, ticker, [0,0,0,0], time_units, p, k, n)  # here you can change the type of mother
        bio_age = result['bio_age']
        age = result['age']
        pad = time_units+1-len(bio_age)  # for avoiding mother matrix to be empty
        bio_age = np.pad(np.array(bio_age, dtype='float'), (0,pad), mode='constant', constant_values=np.nan)
        mothers[trial] = bio_age
        mothers_mean_bio_age = np.nanmean(mothers, axis=0)  # calculating the mean bio age of mother
        mothers_lifespan[trial] = age
    pad = time_units+1-len(mothers_mean_bio_age)
    mean_mothers[iteration] = np.pad(mothers_mean_bio_age, (0, pad), mode='constant', constant_values=np.nan)
    ages_mothers = np.mean(mothers_lifespan)  # chronological mean lifespan of the mothers
total_time = time.process_time() - start


print(f'Population size: {len(population)}')
print('Time for main process:', round(total_time,2), 'seconds.')

print('Average lifespan of mother:', ages_mothers)

df = pd.DataFrame(data=mean_mothers.T, columns=ps)
# print(df)
# print(mothers)

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

ax1 = sn.lineplot(data=df, ax=ax1)  # plotting 3 different iterations
plt.xlabel("Time unit")
plt.ylabel("Mother type")

ax2 = sn.lineplot(data=bio_age, ax=ax2, color='g')  # plotting the single model
plt.xlabel("Time unit")
plt.ylabel("Mother type")

# pprint.pprint(population)
# print(df)
plt.figure(figsize=(8, 6), dpi=300)
# sn.lineplot(data=df)
# sn.heatmap(data=df, cmap='Blues')
plt.show()
print('Finished')