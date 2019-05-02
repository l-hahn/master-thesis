import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import load_iris

def BuildUnsupRfData(XData):
    SynthData = np.append(XData,np.array([ np.random.choice(Row, size=len(Row),replace=True) for Row in np.array(XData).T ]).T,axis=0)
    SynthLabel = np.repeat([0,1], [len(XData),len(XData)],axis=0)
    return SynthData, SynthLabel



#File = open("VarImp.dat","r")
#Data = [list(map(int, Line.rstrip().split(" "))) for Line in File]
#File.close()
#Names = ["A","B","C","D","E"]


#File = open("IrisData","r")
#Data = [list(map(float, Line.rstrip().split(" "))) for Line in File]
#File.close()
#Names = ["Sepal.Length","Sepal.Width","Petal.Length","Petal.Width"]

#SynthData,SynthLabel = BuildUnsupRfData(Data) 
#rnd_clf = RandomForestClassifier(n_estimators = 5000, n_jobs = -1)
#print("Start Cluster")
#rnd_clf.fit(SynthData, SynthLabel)
#print("End Cluster")
#for name, importance in zip(Names, rnd_clf.feature_importances_):
#    print(name, "=", importance)


iris = load_iris()
SynthData,SynthLabel = BuildUnsupRfData(iris.data)
rnd_clf = RandomForestClassifier(n_estimators = 5000, n_jobs = -1)
print("Start Cluster")
rnd_clf.fit(SynthData, SynthLabel)
print("End Cluster")
for name, importance in zip(iris.feature_names, rnd_clf.feature_importances_):
    print(name, "=", importance)



















