import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.metrics import mean_absolute_error


def plot_time_series_class(data, class_name, ax, n_steps=10):
    ''' plot averaged (smoothed out with one standard deviation on top and bottom of it) time series data '''
    time_series_df = pd.DataFrame(data)

    smooth_path = time_series_df.rolling(n_steps).mean()
    path_deviation = 2 * time_series_df.rolling(n_steps).std()

    under_line = (smooth_path - path_deviation)[0]
    over_line = (smooth_path + path_deviation)[0]

    ax.plot(smooth_path, linewidth=2)
    ax.fill_between(
    path_deviation.index,
    under_line,
    over_line,
    alpha=.125)
    
    ax.set_title(class_name)
    

def plot_average_class_representation(dataset, class_dict, label):
    ''' plot averaged (smoothed out with one standard deviation on top and bottom of it) Time Series for each class '''

    # extract class labels
    classes = [v for _, v in class_dict.items()]
    
    # extract class names
    class_names = [k for k, _ in class_dict.items()]
    
    # define subplots
    fig, axs = plt.subplots(
      nrows=len(classes) // 3 + 1,
      ncols=3,
      sharey=True,
      figsize=(14, 8)
    )
    
    for i, cls in enumerate(classes):
        ax = axs.flat[i]
        data = dataset[dataset[label] == cls] \
              .drop(labels=label, axis=1) \
              .mean(axis=0) \
              .to_numpy()
        plot_time_series_class(data, class_names[i], ax)

    fig.delaxes(axs.flat[-1])
    fig.tight_layout()
    
    
class Metrics:
    ''' placeholder of predition metrics '''
    def __init__(self):
        # create separate lists for each metric, test predictions and name
        self.accuracy = []
        self.precision = []
        self.recall =[]
        self.f1 = []
        self.auroc = []
        self.test_pred = []
        self.test_seq = []
        self.names = []
        
    def add_components(self, model_name, X_test_pred, y_test_pred, accuracy, precision, reacall, f1, auroc):
        # add metric components to the appropriate lists
        self.names.append(model_name)
        self.test_pred.append(y_test_pred)
        self.test_seq.append(X_test_pred)
        self.accuracy.append(accuracy)
        self.precision.append(precision)
        self.recall.append(reacall)
        self.f1.append(f1)
        self.auroc.append(auroc)

    def get_metrics(self):
        # return all metric lists
        return self.precision, self.recall, self.f1, self.accuracy, self.auroc
    
    def get_names(self):
        # return list of model names
        return self.names
    
    def get_test_predictions(self):
        # return list of test reconstructed sequences and predictions
        return self.test_seq, self.test_pred
    

def fit_model(model, X_train, epochs, batch_size, verbose, callbacks):
    ''' fit AutoEncoder '''
    # create history of training losses
    history=model.fit(X_train, X_train, epochs=epochs, batch_size=batch_size, verbose=0, callbacks=callbacks)

    # print training history
    plt.figure(figsize=(10,5))
    plt.title('Training history')
    plt.plot(history.history['loss'], label='train loss')
    plt.legend()
    plt.show()
    
    return model    


def make_predictions(model, X_train, X_val, X_test):
    ''' make reconstructions of datasets'''
    
    # make prediction (reconstruction) training set
    X_train_pred = model.predict(X_train)

    # make prediction (reconstruction) validation set
    X_val_pred = model.predict(X_val)

    # make prediction (reconstruction) test set
    X_test_pred = model.predict(X_test)
    
    return X_train_pred, X_val_pred, X_test_pred 


def calculate_reconstrucion_losses(X_train, X_train_pred, X_val, X_val_pred, X_test, X_test_pred):
    ''' calculate reconstruction losses for datasets using mae'''
    
    # calculate reconstrucion loss for each train sequence
    train_loss = [mean_absolute_error( X_train[i].flatten(), X_train_pred[i].flatten() ) for i in range(X_train.shape[0])]

    # calculate reconstrucion loss for each validation sequence
    val_loss = [mean_absolute_error( X_val[i].flatten(), X_val_pred[i].flatten() ) for i in range(X_val.shape[0])]

    # calculate reconstrucion loss for each test sequence
    test_loss = [mean_absolute_error( X_test[i].flatten(), X_test_pred[i].flatten() ) for i in range(X_test.shape[0])]
    
    return train_loss, val_loss, test_loss


def calculate_prediction_metrics(y_test, y_test_pred, verbose=1):
    ''' caclucale and print prediction metrics'''
    
    # calculate test set accuracy prediction
    test_accuracy = metrics.accuracy_score(y_test_pred, y_test)
    test_precision = metrics.precision_score(y_test_pred, y_test)
    test_recall = metrics.recall_score(y_test_pred, y_test)
    test_f1_score = metrics.f1_score(y_test_pred, y_test)
    test_auroc = metrics.roc_auc_score(y_test, y_test_pred)

    if verbose==1:
        print(f'Test precision score: {round(test_precision, 4)}')
        print(f'Test recall score: {round(test_recall, 4)}')
        print(f'Test f1 score: {round(test_f1_score, 4)}')
        print(f'Test accuracy score: {round(test_accuracy, 4)}')
        print(f'Test auroc score: {round(test_auroc, 4)}')
    
    return test_accuracy, test_precision, test_recall, test_f1_score, test_auroc


def select_threshold(train_loss, y_val, val_loss, percentiles):

    # define threshold variable
    best_threshold = None
    best_accuracy = 0

    for percentile in percentiles:
        # set theshold based on percentile on train reconctrucion losses
        testing_threshold = np.percentile(train_loss, percentile)

        # predict validation set classes based on validation reconctrucion losses and testing threshold
        y_val_pred = [int(x < testing_threshold) for x in val_loss]

        # calculate validaion accuracy
        val_accuracy = metrics.accuracy_score(y_val_pred, y_val)

        # compare current validation accuracy with best accuracy
        if val_accuracy > best_accuracy:
            best_accuracy = val_accuracy
            best_threshold = testing_threshold

        print( (f'Perentile:{percentile} Threshold: {round(testing_threshold,4)} ' 
                f'Validation Accuracy: {round(val_accuracy, 4)}' ) )

    print('-'*65)
    print(f'Best validation accuracy: {round(best_accuracy, 4)} for threshold: {round(best_threshold, 4)}')
    
    return best_threshold


def plot_prediction(data, model, title, ax):
    ''' plot simple sequence reconstruction '''
    predictions = model.predict(data)
    pred_losses = mean_absolute_error(predictions.flatten(), data.flatten())
    ax.plot(data.flatten(), label='true')
    ax.plot(predictions.flatten(), label='reconstructed')
    ax.set_title(f'{title} (loss: {np.around(pred_losses, 2)})')
    ax.legend()

    
def plot_recontructions(model, X_1, X_2, threshold, title_1, title_2, nrows=2, ncols=6, figsize=(22, 8) ):
    ''' plot sequences reconstructions '''
    
    fig, axs = plt.subplots( nrows=nrows, ncols=ncols, sharey=True, sharex=True, figsize=figsize)
    fig.suptitle(f'Threshold for anomaly detection: {round(threshold, 4)}\n', fontsize=20, y=1.05)
    
    for i in range(ncols):
        plot_prediction(X_1[i:i+1], model, title=title_1, ax=axs[0, i])

    for i in range(ncols):
        plot_prediction(X_2[i:i+1], model, title=title_2, ax=axs[1, i])

    fig.tight_layout()
    
    
def plot_single_reconstruction(X_true, X_predicted, title, ax):
    ''' plot single sequence reconstruction '''
    pred_losses = mean_absolute_error(X_true.flatten(), X_predicted.flatten())
    ax.plot(X_true.flatten(), label='true')
    ax.plot(X_predicted.flatten(), label='reconstructed')
    ax.set_title(f'{title} (loss: {np.around(pred_losses, 2)})')
    ax.legend()


def compare_reconstructions(X_test, y_test, model_names, sequences, predictions, title, target1, target2, figsize):
    ''' compare true sequences with reconstructed sequences for all models'''
    
    # caclulate nrows based on number of models passed to the function
    nrows = len(model_names) 
    
    fig, axs = plt.subplots( nrows=nrows, ncols=6, sharey=True, sharex=True, figsize=figsize )
    fig.suptitle(t=title, fontsize=20, y=1.02)
    
    for row in range(nrows):
        # extract predictions for particular model
        y_pred = np.array(predictions[row])
        
        # create subests to compare
        X_pred = sequences[row][(y_test==target1) & (y_pred==target2)]
        X_true = X_test[(y_test==target1) & (y_pred==target2)]
        
        # set number of cols no more than 6
        ncols = min(len(X_true), 6)
        
        # compare single true sequence with reconstruction
        for col in range(ncols):
            plot_single_reconstruction(X_true[col:col+1], X_pred[col:col+1], title=model_names[row], ax=axs[row, col])    

    fig.tight_layout()


def plot_results(results):
    results.plot.barh(figsize=(13,5), xlim=(0.93, 0.99))
    plt.axvline(x=0.94, color='k', alpha=1, ls='--')
    plt.axvline(x=0.95, color='k', alpha=1, ls='--')
    plt.axvline(x=0.96, color='k', alpha=1, ls='--')
    plt.axvline(x=0.97, color='k', alpha=1, ls='--')
    plt.axvline(x=0.98, color='k', alpha=1, ls='--')
    plt.legend(loc='lower right')
    plt.title('Comparission of all results')
    plt.tight_layout()
    plt.show()