import keras.backend as K
import re


def total_accuracy(y_true, y_pred):
    return K.sum(K.cast(K.equal(K.argmax(y_true, axis=-1), K.argmax(y_pred, axis=-1)), 'float32') * K.sum(y_true,
                                                                                                          axis=-1)) / K.sum(
        y_true)


INTERESTING_CLASS_ID = 1  # Choose the class of interest


def sensitivity(y_true, y_pred):
    tr = K.argmax(y_true, axis=-1)
    te = K.argmax(y_pred, axis=-1)
    lp = K.cast(K.sum(y_true, -1), 'float32')
    sc = K.cast(K.equal(tr, 1), 'float32') * lp
    sc1 = K.cast(K.equal(te, 1), 'float32') * lp
    tot = K.cast(K.sum(sc), 'float32')
    lp = K.cast(K.sum(y_true, -1), 'float32')
    return K.sum(sc * sc1) / (tot + K.epsilon())


def precision(y_true,y_pred):
    tr = K.argmax(y_true, axis=-1)
    te = K.argmax(y_pred, axis=-1)
    lp = K.cast(K.sum(y_true, -1), 'float32')
    sc = K.cast(K.equal(tr, 1), 'float32')*lp
    sc1 = K.cast(K.equal(te, 1), 'float32')*lp
    tot = K.cast(K.sum(sc1), 'float32')
    lp = K.cast(K.sum(y_true,-1), 'float32')
    return K.sum(sc*sc1)/(tot + K.epsilon())


def f1(y_true,y_pred):
    tr = K.argmax(y_true, axis=-1)
    te = K.argmax(y_pred, axis=-1)
    lp = K.cast(K.sum(y_true, -1), 'float32')
    sc = K.cast(K.equal(tr, 1), 'float32')*lp
    sc1 = K.cast(K.equal(te, 1), 'float32')*lp
    tot = K.cast(K.sum(sc1), 'float32')
    lp = K.cast(K.sum(y_true, -1), 'float32')
    tot1 = K.cast(K.sum(sc), 'float32')
    return 2/(1/(K.sum(sc*sc1)/(tot + K.epsilon()))+1/(K.sum(sc*sc1)/(tot1 + K.epsilon())))


def calc_sov(true, pred, label='1', calc_sigma=True):
    """
    Calculates segment overlap based on Zemla et al. Proteins 1999, 34:220–223
    :param true: true assignment
    :param pred: predicted assignment
    :param label: label for which SOV is calculated
    :param calc_sigma: determines whether sigma factor is calculated
    :return SOV value for prootein:
    """
    assert len(true) == len(pred)
    if true.count(label) > 0:
        Ncc = 0
        true_segments = []
        c = 0
        for match in re.finditer(r"([%s]+)" % label, true):
            true_segments.append([])
            for i in range(match.start(), match.end()):
                true_segments[c].append(i)
            c += 1
        c = 0
        pred_segments = []
        for match in re.finditer(r"([%s]+)" % label, pred):
            pred_segments.append([])
            for i in range(match.start(), match.end()):
                pred_segments[c].append(i)
            c += 1
        summ = 0
        for segment in true_segments:
            matching_segment = False
            for segment2 in pred_segments:
                minov = len((set(segment) & set(segment2)))
                maxov = len((set(segment) | set(segment2)))
                if minov > 0:
                    Ncc += len(segment)
                    matching_segment = True
                    len_s1_adj = int(0.5 * len(segment))
                    len_s2_adj = int(0.5 * len(segment2))
                    sigma = min(maxov - minov, minov, len_s1_adj, len_s2_adj)
                    if not calc_sigma:
                        sigma = 0
                    summ += float(len(segment) * (minov + sigma) / maxov)
            if not matching_segment:
                Ncc += len(segment)
        return (float(summ / Ncc))

