import matplotlib.pyplot as plt
import numpy as np
import sys
import json

LOG_RATE = 20

def plot(runcard, label):
    rep_list = []
    # open every replica file and add the chi^2exps.log dicts to a list
    all_chi2 = []
    time = {'walltime':[], 'cputime':[]}
    for replica in range(1,130):
        text = ''
        try:
            with open(f'{runcard}/nnfit/replica_{replica}/chi2exps.log', 'r') as f:
                for line in f:
                    text = text + line.strip()
                rep_list.append(eval(text))
                
            with open(f'{runcard}/nnfit/replica_{replica}/testchi^2.txt', 'r') as f:
                for line in f:
                    all_chi2.append(float(line.strip()))
                    
            with open(f'{runcard}/nnfit/replica_{replica}/{runcard}.json') as f:
                data = json.load(f)
                time['walltime'].append(data['timing']['walltime']['Total'])
                time['cputime'].append(data['timing']['cputime']['Total'])

        except FileNotFoundError:
            break
        except:
            print('Something else went wrong')

    # make 2 arrays with shape (amount of reps, # of epochs/log_rate)
    longest_rep = len(max(rep_list, key=len))
    training = np.zeros((len(rep_list), longest_rep))
    validation = np.zeros((len(rep_list), longest_rep))
    times = [sum(time[t])/len(time[t]) for t in time]

    # add data to both arrays from log file
    epochs = []
    for i,rep in enumerate(rep_list):
        last_epoch = 0
        for j, epoch in enumerate(rep):
            training[i, j] = rep[epoch]['total']['training']
            validation[i, j] = rep[epoch]['total']['validation']
            # if runcard == 'runcard_kernels':
            #     if float(training[i,j]) < 1.05 and int(epoch) ==5000:
                    # print(i, training[i,j])
            last_epoch = epoch
        epochs.append(int(last_epoch))
        
        #if rep has less epochs than longest_rep, fill those epoch values with the final loss value
        if j < longest_rep-1:
            training[i,j:] = training[i,j]
            validation[i,j:] = validation[i,j]
    epoch_mean =  sum(epochs)/len(epochs)

    # find mean and 95% loss interval for the training and validation error
    train_avg = np.mean(training, axis=0)
    train_ci = 1.96 * np.std(training, axis=0)/np.sqrt(training.shape[0])
    val_avg = np.mean(validation, axis=0)
    val_ci = 1.96 * np.std(validation, axis=0)/np.sqrt(validation.shape[0])

    x = np.arange(LOG_RATE, LOG_RATE * longest_rep + LOG_RATE, LOG_RATE)

    colors = {
        'runcard_no_mod':'darkslategrey',
        'runcard_skip_connections':'indianred',
        'runcard_skip_full': 'sandybrown',
        'runcard_kernels': 'cadetblue'
    }
    plt.plot(x, train_avg, label = f'{label}', color= colors[runcard])
    # plt.fill_between(x, train_avg - train_ci, train_avg + train_ci, alpha=0.5, label = f'{label}')
    # plt.plot(x, val_avg, label = f'{runcard} Total validation')
    # plt.fill_between(x, val_avg - val_ci, val_avg + val_ci, alpha=0.5)
    plt.ylim(1,1.5)
    plt.legend()
    plt.xlabel('Epoch')
    plt.ylabel(r'$\chi^2$')
    return all_chi2, epoch_mean, times

    
if __name__ == "__main__":
    chi=[]
    means=[]
    times =[]
    runcards= ['runcard_no_mod', 'runcard_skip_connections', 'runcard_skip_full', 'runcard_kernels']
    labels = ['No Modifications', 'Skip Connections', 'Fully Connected', 'Kernel Methods']
    for runcard, label in zip(runcards, labels):
        res = plot(runcard, label)
        chi.append(res[0])
        means.append(res[1])
        times.append(res[2])
    plt.title('Average training loss')
    plt.savefig(f'plot_losses.png')
    
    plt.figure()
    plt.ylim(1.05,1.15)
    plt.boxplot(chi,labels=labels, 
                patch_artist = True, 
                boxprops = dict(facecolor = "darkslategrey", alpha=0.5), 
                medianprops = dict(color = "indianred", linewidth = 1.5))
    plt.xticks(rotation=25, ha='right')
    plt.ylabel(r'$\chi^2$')
    plt.title('Average test loss')
    plt.tight_layout()
    plt.savefig(f'plot_chiexp.png')
    
    fig, ax = plt.subplots()
    plt.title('Epochs and timing')
    ax2 = ax.twinx()
    xticks = range(len(labels))
    ax.bar([x -0.2 for x in xticks], means, width=0.2, color = 'darkslategrey')
    ax2.bar(xticks, [t[0] for t in times], width=0.2, color='indianred')
    # plt.legend(loc='upper left')

    # giving labels to the axises
    ax.set_ylabel('# of epochs', color = 'darkslategrey')
    ax.set_ylim(4000,5000)
    
    # secondary y-axis label
    ax2.set_ylabel('Wall Time (s)', color = 'indianred')
    plt.tight_layout()
    ax.set_xticks(xticks, labels, rotation=25, ha='right')
    ax.tick_params(axis='y', labelcolor='darkslategrey')
    ax2.tick_params(axis='y', labelcolor='indianred')
    plt.tight_layout()
    plt.savefig(f'plot_epochs.png')
    
