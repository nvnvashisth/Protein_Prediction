import gzip
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report


def generate_fasta_reader(filename, gzipped=False, strip=True):
    """This module function is used to create a generator function to read various flavors
    of fasta formated data"""

    def getFastaEntry():
        filehandle = None
        if gzipped == True:
            filehandle = gzip.open(filename, 'rb')
        else:
            filehandle = open(filename, 'rb')

        header = ""
        body = ""
        fileStart = True
        for line in filehandle:
            # convert from binary to string
            line = line.decode()
            # remove trailing white space if strip == True
            if strip:
                line = line.strip()
            # line is comment
            if line.startswith('#'):
                continue
            # line is empty
            if len(line) == 0:
                continue
            # line is first header
            if (fileStart == True and line.startswith('>')):
                header = line
                fileStart = False
                continue
            # line is header but not first
            if line.startswith('>'):
                ret = (header, body)
                header = line
                body = ""
                # print("ejecting fasta entry")
                yield ret
            else:
                body = body + line

        # last entry / cleanup
        if fileStart == False:
            filehandle.close()
            yield (header, body)

    return getFastaEntry


class Feature_Generator:
    """
    The class Feature_Generator uses the method to create a fasta reader function to read in the data. The layout of
    the data is header line, sequence line, label line. It requires a file name for the input and half of the window size.
    This is to ensure a symmetric window with the residue to describe in the center and avoid problems with even numbers
    as window sizes. The __init__ function takes (besides self) two parameters filename and half_window_size. The final
    window size is 2 * half_window_size + 1. The function __sliding_window is called by __init__
    :filename: the input filename
    :half_window_size: the number of neighboring residues on each side of the central residue. Default to 3.
    """

    def __init__(self, filename, half_window_size=3):
        """
        Needs an input file name and half_window_size
        """
        #TODO
        self.half_window_size = half_window_size
        self.fasta_list = self.fasta_reader(filename)
        self.all_residues = list()
        self.__sliding_window(self.fasta_list,half_window_size)
        self.df = None

    def get_fasta_entry(self,filename, gzipped=False, strip=True):
        filehandle = None
        if gzipped == True:
            filehandle = gzip.open(filename, 'rb')
        else:
            filehandle = open(filename, 'rb')

        header = ""
        body = ""
        fileStart = True
        for line in filehandle:
            # convert from binary to string
            line = line.decode()
            # remove trailing white space if strip == True
            if strip:
                line = line.strip()
            # line is comment
            if line.startswith('#'):
                continue
            # line is empty
            if len(line) == 0:
                continue
            # line is first header
            if (fileStart == True and line.startswith('>')):
                header = line
                fileStart = False
                continue
            # line is header but not first
            if line.startswith('>'):
                ret = (header, body)
                header = line
                body = ""
                # print("ejecting fasta entry")
                yield ret
            else:
                body = body + line

        # last entry / cleanup
        if fileStart == False:
            filehandle.close()
            yield (header, body)


    def fasta_reader(self,filename):

        fasta_list = list()
        try:
            fasta_func = self.get_fasta_entry(filename)

            #valid = True
            while(True):
                fasta_val = fasta_func.__next__()
                length = int(len(fasta_val[1])/2)
                protein_seq = fasta_val[1][0:length]
                label_seq = fasta_val[1][length:]
                fasta_list.append((fasta_val[0],protein_seq,label_seq))

        except Exception as e:
            print("All read")

        return fasta_list



    def get_feature_table(self):
        """
        Return an unmodified version of the encoding as list of lists
        :return: Returns self.all_residues
        """
        return self.all_residues

    def __sliding_window(self, fasta_list, half_window_size=3):
        """
        The function sliding_window expects a protein sequence string and a residues' label string of the same length.
        The parameter half_window_size determines how many residues on each side of a central residue are included.
        Sequence residues at [0,window_size]  or [len(sequence)-window_size, len(sequence)] are not recoded and omitted
        in the result value. For every eligible residue in the sequence there will be a list consisting of the central
        residue, the neighbors (from half_window_size left of the central residue to half_window_size right of the
        central residue), and the state label. Such a list is added to self.residues for every recoded residue, i.e. upon
        completion self.all_residues is a list of lists.
        """
        #TODO
        for fasta in fasta_list:
            seq = fasta[1]
            label = fasta[2]
            for i in range(half_window_size,len(seq)-half_window_size):
                residue_attr = list()
                # Residue itself
                residue_attr.append(seq[i])

                # For preceding residues
                for j in range(i-half_window_size,i):
                    residue_attr.append(seq[j])

                # For post sequences
                for j in range(i+1,i+half_window_size+1):
                    residue_attr.append(seq[j])

                # For the label
                residue_attr.append(label[i])
                self.all_residues.append(residue_attr)

    def create_dataframe(self):
        """
        Creates a pandas.DataFrame from the data stored in self.all_residues for easier processing and stores it in self.df.
        It adds meaningful column names so that the naming scheme is:
        center, -half_win_size, ..., -1, +1, +half_win_size, state

        This DataFrame is stored for further processing
        """
        self.df = pd.DataFrame.from_records(self.all_residues)
        # add code to give meaningful columns names, including the one base on win_size, here
        # TODO
        window_size = self.half_window_size
        new_columns = ["center"]
        # For negative values
        neg_val = -1*window_size
        for i in range(neg_val,0):
            new_columns.append(str(i))

        # For positive values
        for i in range(1,window_size+1):
            new_columns.append(str(i))

        new_columns.append("state")
        self.df.columns = new_columns

    def get_data_frame(self):
        """
        Returns the data as DataFrame
        :return: self.df
        """
        if self.df is None:
           self.create_dataframe()
        return self.df

    def remove_rows_from_df(self, col, val):
        """
        Remove row from the DataFrame df, where columns col has value val
        :param col: column which is test for value val
        :param val: filtered value
        :return: No return, modifies self.df, only lines where df["col"] != val are kept
        """
        #TODO
        index = 0
        for i,column in enumerate(self.df.columns):
            if col == column:
                index = i
                break

        val = self.df[self.df.columns[index]] != val
        self.df = self.df[val]

    def modify_cells_from_df(self, col, old_val, new_val):
        """
        Change the value of cells if DataFrame df specified by column col and value old_val to value new_val
        :param col: column where to look for old_val
        :param old_val: value to be replaced
        :param new_val: replacement value
        :return: No return, modifies self.df in place
        """
        #TODO
        index = 0
        for i, column in enumerate(self.df.columns):
            if col == column:
                index = i
                break
        self.df[self.df.columns[index]].replace(old_val,new_val, inplace=True)



    def create_full_binarized_dataframe(self):
        """
        This function returns a version of a pandas.DataFrame where all nominal attributes have
        been converted into binary ones  including the class label
        :return: pandas.DataFrame
        """
        #TODO
        df_full_bin = pd.get_dummies(self.df)
        return df_full_bin

    def create_partial_binarized_dataframe(self, list_of_features):
        """
        This function returns a version of a pandas.DataFrame where all specified nominal
        attributes have been binarized
        """

        df_part_bin = pd.get_dummies(self.df, columns=list_of_features)
        return df_part_bin

    def encode_multi_class_label(self, df_in, class_label):
        """
        This function modifies the given data frame by encoding a nominal class attribute.
        The class labels are converted to values between 0 and n-1 in place if you have n classes.
        :param: this is the input pandas.DataFrame which class attribute you want to modify
        :param class_label: name of the class attribute
        """
        le = LabelEncoder()
        le.fit(df_in[class_label])
        # print(le.classes_)
        df_in[class_label] = le.transform(df_in[class_label])




class Machine_Learner:
    """
    This class creates and use a LogisticRegression classifier
    """
    def __init__(self, seed):
        self.seed = seed
        self.clf = LogisticRegression(random_state=seed)


    def get_clf(self):
        return self.clf

    def set_df(self, df, class_label):
        self.df = df
        self.class_label = class_label

    def initialize_data(self):
        """
        Scikit machine learning algorithms typically take the feature data a 2D numpy array and the label data as 1D numpy array.
        :return: no return value
        """
        data = self.df.as_matrix()
        rows = data.shape[0]
        cols = data.shape[1]
        self.X = data[:,1:]
        self.y = data[:,0]
        self.labels = np.unique(self.y)

    def get_X(self):
        return self.X

    def get_y(self):
        return self.y

    def initialize_training_test_sets(self):
        #TODO
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(self.X, self.y,random_state=self.seed)

    def get_training_sets(self):
        return(self.X_train, self.y_train)

    def get_test_sets(self):
        return(self.X_test, self.y_test)

    def train_classifier(self):
        #TODO
        self.clf.fit(self.X_train, self.y_train)

    def make_predictions(self):
        #TODO
        return self.clf.predict(self.X_test)

    def get_confusion_matrix(self):
        return confusion_matrix(self.y_test, self.make_predictions(),self.labels)

    def get_classification_report(self):
        return classification_report(self.y_test, self.make_predictions(),self.labels)



from collections import Counter

def main():
    fgen = Feature_Generator("data/opm_tmps_stud.fasta", half_window_size=3 )
    fgen.create_dataframe()
    fgen.remove_rows_from_df("state", "U")
    fgen.modify_cells_from_df("state", 'h', 'H')
    cols = fgen.get_data_frame().columns
    cols = cols[:-1]
    df = fgen.create_partial_binarized_dataframe(list(cols))
    fgen.encode_multi_class_label(df, "state")

    clf = Machine_Learner(42)
    clf.set_df(df, "state")
    clf.initialize_data()
    clf.initialize_training_test_sets()
    clf.train_classifier()
    print(clf.get_confusion_matrix())
    print(clf.get_classification_report())

if __name__ == "__main__":
    # execute only if run as a script
    #
    main()
