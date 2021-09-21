'''
Spam predictor

pip install -U scikit-learn

X are the data features

y are the labels [ 0, 0, 1 ...]

Download the enron database for testing

http://www2.aueb.gr/users/ion/data/enron-spam/

Should work on any of the datasets:

seq 1 6 | parallel -j 1 wget -q -nc http://www.aueb.gr/users/ion/data/enron-spam/preprocessed/enron{}.tar.gz

'''
import logging
import sys, os

import plac
from joblib import dump, load

logger = logging.getLogger("engine")

try:
    from sklearn.feature_extraction.text import CountVectorizer
    from sklearn.metrics import classification_report
    from sklearn.model_selection import train_test_split
    from sklearn.naive_bayes import MultinomialNB
    from sklearn.pipeline import make_pipeline
    has_sklearn = True
except ImportError as exc:
    logger.error("sklearn not installed, no predictions are generated")
    has_sklearn = False

def load_model(model="spam.model"):
    nb = load(model)
    return nb


def classify_content(content, model):
    """
    Classify content
    """
    if not has_sklearn:
        return 0

    try:
        nb = load_model(model)
        y_pred = nb.predict([content])
    except Exception as exc:
        logger.error(exc)
        y_pred = [0]

    return y_pred[0]


def fit_model(X, y):

    nb = make_pipeline(

        CountVectorizer(),

        MultinomialNB(),

        # LinearSVC(),

        # RandomForestClassifier(),

    )

    nb.fit(X, y)

    return nb



def evaluate_model(fname, model):

    X, y = parse_file(fname=fname)

    X_train, X_test, y_train, y_test = train_test_split(X, y)

    if model:
        nb = load_model(model)
    else:
        nb = fit_model(X_train, y_train)

    y_pred = nb.predict(X_test)
    rep = classification_report(y_test, y_pred)

    print(rep)


def parse_file(fname):
    import tarfile

    # Take filename from command line
    tar = tarfile.open(name=fname, mode='r:gz', fileobj=None)
    elems = filter(lambda t: t.isreg(), tar)
    X, y = [], []
    for info in elems:
        # Fill in the labels
        y.append(int("spam" in info.name))
        stream = tar.extractfile(info)
        content = stream.read().decode("utf-8", errors="ignore")
        X.append(content)

    spam_count = sum(filter(None, y))
    ham_count = len(y) - spam_count

    logger.info(f"parsed: {ham_count} ham, {spam_count} spam")

    return X, y


def build_model(fname, model):
    '''
    wget -nc http://www.aueb.gr/users/ion/data/enron-spam/preprocessed/enron1.tar.gz
    '''

    # Generate features.
    X, y = parse_file(fname)

    # Fits the model.
    nb = fit_model(X, y)

    logger.info(f"fitted model to: {fname}")

    # Save the model.
    if model:
        logger.info(f"saving model to: {model}")
        dump(nb, model)

    return nb

    # evaluate_model(X, y, model=model)

@plac.pos('fname')
@plac.flg('build')
@plac.flg('eval_', help="evaluate model ")
@plac.opt('model')
@plac.flg('classify')
def main(classify, build, model, eval_, fname):

    if build:
        build_model(fname=fname, model=model)

    if eval_:
        evaluate_model(fname=fname, model=model)

    if classify:
        content = open(fname, 'rt').read()
        res = classify_content(content=content, model=model)
        print ("spam" if res else "ham")

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)

    plac.call(main)
