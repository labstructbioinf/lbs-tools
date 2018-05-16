from keras import callbacks
from sklearn.metrics import roc_curve, auc, precision_score, recall_score, f1_score
import numpy as np


class Logger(callbacks.Callback):
    def __init__(self, lab_pos=2,out_path='', out_fn='best_model.h5', out_log = ''):
        super(Logger, self).__init__()
        self.lab_pos = lab_pos
        self.f1 = 0
        self.auc = 0
        self.path = out_path
        self.fn = out_fn
        self.log = out_log

    def on_train_begin(self, logs={}):
        return

    def on_train_end(self, logs={}):
        f = open(self.log, 'a')
        f.write('%s\n' % (self.f1))
        f.close()

    def on_epoch_begin(self, epoch, logs={}):
        return

    def on_batch_begin(self, batch, logs={}):
        return

    def on_batch_end(self, batch, logs={}):
        return

    def on_epoch_end(self, epoch, logs={}):
        cv_pred = self.model.predict(self.validation_data[0])
        scr_val, wyb_val, tr_val = [], [], []
        for j, j1 in zip(self.validation_data[1], cv_pred):
            for i, i1 in zip(j, j1):
                if max(i) == 1:
                    tr_val.append(np.argmax(i))
                    wyb_val.append(np.argmax(i1))
                    scr_val.append(i1[self.lab_pos])
        rec_val = recall_score(tr_val, wyb_val, average='micro', labels=[self.lab_pos])
        prec_val = precision_score(tr_val, wyb_val, average='micro', labels=[self.lab_pos])
        f1_val = f1_score(tr_val, wyb_val, average='micro', labels=[self.lab_pos])
        #fpr, tpr, thresholds = roc_curve(tr_val, scr_val, pos_label=self.lab_pos)
        #auc_val = auc(fpr, tpr)
        #print('Val F1_score ' + '%.3f' % f1_val + ' ' + 'Val sens: ' + '%.3f' % rec_val + ' Val prec: ' + '%.3f' % prec_val + ' AUC: ' + '%.3f' % auc_val)
        if self.f1 < f1_val:
            f1_val2 = '%.3f' % f1_val
            rec_val2 = '%.3f' % rec_val
            prec_val2 = '%.3f' % prec_val
            print("Best F1 score: %s (prec: %s, sens: %s)" % (f1_val2, prec_val2, rec_val2))
            self.f1 = f1_val
            self.model.save(self.path + self.fn, overwrite=True)
            #self.model.save(self.path + "f1_" + str(epoch) + '_' + str(rec_val2) + '_' + str(prec_val2) + '_' + str(
                    #f1_val2) + ".h5", overwrite=True)
        #if self.auc < auc_val:
            #pass
            #auc_val2 = '%.3f' % auc_val
            #print("Best AUC score!")
            #self.auc = auc_val
            #self.model.save(self.path + "auc_" + str(epoch) + '_' + str(auc_val2)  + ".h5", overwrite=True)
        return
