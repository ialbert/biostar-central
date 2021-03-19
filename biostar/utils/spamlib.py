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

from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from sklearn.pipeline import make_pipeline

from sklearn.feature_extraction.text import CountVectorizer
from sklearn.feature_extraction.text import TfidfTransformer

from sklearn.naive_bayes import MultinomialNB, GaussianNB, BernoulliNB
from sklearn.svm import SVC, NuSVC, LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

from sklearn.base import TransformerMixin
from itertools import *
import sys

from joblib import dump, load


def load_model(fname="spam.model"):
    nb = load(fname)
    return nb


def build_model(X, y, fname="spam.model"):
    nb = make_pipeline(
        CountVectorizer(),
        MultinomialNB(),
        #LinearSVC(),
        #RandomForestClassifier(),
    )
    nb.fit(X, y)
    if fname:
        dump(nb, fname)
    return nb


def evaluate_model(X, y, model=None):
    X_train, X_test, y_train, y_test = train_test_split(X, y)
    if model:
        nb = load_model(model)
    else:
        nb = build_model(X_train, y_train)
        dump(nb, "spam.model")

    y_pred = nb.predict(X_test)
    rep = classification_report(y_test, y_pred)
    print(rep)

def parse_file(stream):
    for line in stream:
        for word in line.split():
            word = word.decode("utf-8", errors="ignore").strip()
            yield word


def test(fname=None, model=None):
    '''
    wget -nc http://www.aueb.gr/users/ion/data/enron-spam/preprocessed/enron1.tar.gz
    '''
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
    print(f"Posts: {ham_count} ham, {spam_count} spam")

    evaluate_model(X, y, model=model)


def main():
    fname = sys.argv[1]

    if len(sys.argv) > 2:
        model = sys.argv[2]
    else:
        model = None

    test(fname=fname, model=model)



if __name__ == '__main__':
    main()
