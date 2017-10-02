
# coding: utf-8

# In[1]:




# In[14]:

import csv

class Entry:
    def __init__(self, protein_db_id="", protein_acc="", gene_db_id="", gene_name="", sec_id="", synonyms = []):
        self.gene_db_id = gene_db_id
        self.gene_name = gene_name
        self.protein_acc = protein_acc
        self.protein_db_id = protein_db_id
        self.sec_id = sec_id
        self.synonyms= synonyms
    def set_synonyms(self, synonyms):
        self.synonyms = synonyms
        
class Log:
    def __init__(self, replacement_string=None, status=None, message_id=None, message=""):
        self.replacement_string = replacement_string
        self.status = status
        self.message_id = message_id
        self.message = message
        
class Preferences:
    # CHANGE
    def __init__(self, prefer_small_gene_names=None, prefer_first_selection_on_multiple=None, prefer_remember_selection=None, prefer_output=None, prefer_gene_db_id=False):
        self.prefer_small_gene_names = prefer_small_gene_names 
        self.prefer_first_selection_on_multiple = prefer_first_selection_on_multiple
        self.prefer_remember_selection = prefer_remember_selection
        self.prefer_output = prefer_output
        self.prefer_gene_db_id = prefer_gene_db_id
        
class globalMessageId():
    single_replacement = "SINGLE_REPLACEMENT"
    single_shortest = "SINGLE_SHORTEST"
    single_longest = "SINGLE_LONGEST"
    single_cg = "SINGLE_CG"
    multiple_replacement = "MULTIPLE_REPLACEMENT"
    multiple_shortest = "MULTIPLE_SHORTEST"
    multiple_longest = "MULTIPLE_LONGEST"
    multiple_cg = "MULTIPLE_CG"
    no_result = "NO_RESULT"
    
class globalSelectionHistory():
    history_list = []

def set_log(status, replacement_string, message_id, message):
    return Log(replacement_string, status, message_id, message)

gene_container = list()

def read_id_file(filenames):
    for file in filenames:
        with open(file, 'rt') as tsvin:
            tsvin = csv.reader(tsvin, delimiter='\t')    
            current_sec_id = ""
            currentEntry = None
            for row in tsvin:
                sec_id = ""
                try:
                    sec_id = row[4]
                except IndexError:
                    pass
                
                if current_sec_id != sec_id: # We are on new entry
                    #Lets add the last entry, if its not the first
                    if currentEntry != None:
                        gene_container.append(currentEntry)
                    
                    #Lets add initial data
                    currentEntry = Entry(row[0], row[1], row[2], row[3], current_sec_id)
                    synonyms = []
                    synonyms.append(row[5])
                    currentEntry.set_synonyms(synonyms)
                    current_sec_id = sec_id
                else: # We are on new synonym
                    currentEntry.synonyms.append(row[5])
                
                
                
def find_gene_name(term, filenames=[], debug=True):
    if len(gene_container) == 0:
        if debug: print("INFO: Loading " + str(filenames))
        read_id_file(filenames)
        
    # CHANGE
    matches = [(entry.gene_name, entry.sec_id, entry.gene_db_id) for entry in gene_container if entry_exist(entry, term)]
    if len(matches) == 0:
        result = ([],[])
    else:
        gene_names, sec_ids, gene_db_ids = zip(*matches)
        #CHANGE
        result= (list(set([gene_name for gene_name in gene_names if gene_name.strip() != ""])), list(set([sec_id for sec_id in sec_ids if sec_id.strip() != ""])), list(set([gene_db_id for gene_db_id in gene_db_ids if gene_db_id.strip() != ""])))    
    
    if len(result[0]) == 0:
        if debug: print('WARNING: No matches found for term="' + term + '"')
    elif debug and len(result[0]) > 1:
        if debug: print('WARNING: ' + str(len(result[0])) + '  matches found for term="' + term + '"')
    
    return result
    
def get_user_input(input_term, results, preferences):
    if preferences.prefer_first_selection_on_multiple: # if user preferred first result
        return results[0] # return first result
    if preferences.prefer_remember_selection: # if user wanted to remember their selections
        for history_results in globalSelectionHistory.history_list: # see if the input results were seen before
            if results == history_results[0]: # if they were
                return history_results[1] # return the previous selection
    
    # Print to client
    
    print("More than 1 legitimate values found for: " + str(input_term))
    print("Options: ")
    for i in range(0, len(results)):
        print(str(i) + ": " + results[i])
        
    while True:
        user_input = input("Type index of desired result, or type custom name: ") # get user selection
        input_index = None
        selection = None
        try:
            input_index = int(user_input) # Try to convert to an integer
            if input_index in range(0, len(results)): # if the input converted successful, see if the integer is valid
                selection = results[input_index] # if it is valid, than set the selection as the proper input
            else: # if it is not valid (out of range)
                print("Must be an index shown above. Please try again.") # tell user to try again
        except (TypeError, ValueError) as err:  # if the input is not an integer
            selection =  user_input # set the string as the selection
        if preferences.prefer_remember_selection: # if user wanted to remember their selections
            globalSelectionHistory.history_list.append((results, selection)) # append to the history: results and selection
        if selection is not None: # if selection was made
            return selection # return the selection
        # else redo the prompt

        
def retrieve_gene_name_facilitator(input_term, input_filenames, preferences):
    try:
        results = find_gene_name(term=input_term, filenames=input_filenames, debug=preferences.prefer_output) # Get results from converter
        log = Log()
        if preferences.prefer_gene_db_id == True:
            if len(results[2]) == 1:
                replacement_string = results[2][0]
                log = set_log("info", replacement_string, globalMessageId.single_replacement, str(input_term) + " -> " + replacement_string) # create the log
            if len(results[2]) > 1:
                replacement_string = get_user_input(input_term, results[2], preferences) # prompt the user with all gene names and set their slection as the replacement string
                log = set_log("warning", replacement_string, globalMessageId.multiple_replacement, str(input_term) + " -> " + replacement_string) # create the log
        elif len(results[0]) == 1: # if there is 1 gene name given
            replacement_string = results[0][0] # set it as the replacement string
            log = set_log("info", replacement_string, globalMessageId.single_replacement, str(input_term) + " -> " + replacement_string) # create the log
        elif len(results[0]) > 1: # if there is more than 1 gene name given
            if preferences.prefer_small_gene_names == True: # if the user decided to take the shortest gene name
                new_results = [name for name in results[0] if len(name) == len(min(results[0], key=len))] # gets all shortest gene names (must be at least 1)
                if len(new_results) == 1: # if there is only one shortest gene name
                    replacement_string = new_results[0] # set it as the replacement string
                    log = set_log("info", replacement_string, globalMessageId.single_shortest, str(input_term) + " -> " + replacement_string) # create the log
                else: # if there is more than 1
                    replacement_string = get_user_input(input_term, new_results, preferences) # prompt the user and set their selection as the replacement string
                    log = set_log("warning", replacement_string, globalMessageId.multiple_shortest, str(input_term) + " -> " + replacement_string) # create the log                        
            elif preferences.prefer_small_gene_names == False: # if the user decided to take the longest gene name
                new_results = [name for name in results[0] if len(name) == len(max(results[0], key=len))] # gets all longest gene names (must be at least 1)
                if len(new_results) == 1: # if there is only one longest gene name
                    replacement_string = new_results[0] # set it as the replacement string
                    log = set_log("info", replacement_string, globalMessageId.single_longest, str(input_term) + " -> " + replacement_string) # create the log
                else: # if there is more than 1
                    replacement_string = get_user_input(input_term, new_results, preferences) # prompt the user and set their selection as the replacement string
                    log = set_log("warning", replacement_string, globalMessageId.multiple_longest, str(input_term) + " -> " + replacement_string) # create the log
            else: # if the user decided to manual select the multiples
                replacement_string = get_user_input(input_term, results[0], preferences) # prompt the user with all gene names and set their slection as the replacement string
                log = set_log("warning", replacement_string, globalMessageId.multiple_replacement, str(input_term) + " -> " + replacement_string) # create the log
        elif len(results[0]) == 0: # if there is no gene names
            if len(results[1]) == 1: # but there is 1 CG id
                replacement_string = results[1][0] # set it as the replacement
                log = set_log("warning", replacement_string, globalMessageId.single_cg, str(input_term) + " -> " + replacement_string) # create the log
            elif len(results[1]) > 1: # if, for some reason, there are more than 1 CG id
                print("No gene names found, but multiple valid CG strings found - follow the instructions below:") # let the user know that the options are cg's
                replacement_string = get_user_input(input_term, results[1], preferences) # prompt the user and set their selection as the replacement string
                log = set_log("warning", replacement_string, globalMessageId.multiple_cg,str(input_term) + " -> " + replacement_string) # create the log
            elif len(results[1]) == 0: # if there are also no CG ids
                replacement_string = str(input_term) # Then cannot replace, just set the input FBgn as the replacement string
                log = set_log("critical", replacement_string, globalMessageId.no_result, str(input_term) + " -> " + replacement_string) # create the log
        
        row = [log.message_id, log.message, str(results[0]), str(results[1])] # set the row to append
        if log.status is "info" or log.status is "critical" or log.status is "warning":
            if preferences.prefer_output: print(row) # print row if critical or warning
        return log.replacement_string # return the replacement string
    except IndexError: 
        return input_term

def entry_exist(entry, term):
        return entry.gene_db_id.upper() == term.upper() or entry.protein_acc.upper() == term.upper() or entry.protein_db_id.upper() == term.upper() or entry.gene_name.upper() == term.upper() or term.upper() in [syn.upper() for syn in entry.synonyms]




# In[2]:

# TEST
if __name__ == "__main__":
    #print(retrieve_gene_name_facilitator(input_term="FBgn0053855", input_filenames=["../../res/mortimer_gene_ids.txt", "../../res/flymine_id_list_3.tsv"], 
    #                    preferences=Preferences(prefer_small_gene_names=True, prefer_first_selection_on_multiple=True, prefer_remember_selection=True, prefer_output=False)))
    print(retrieve_gene_name_facilitator(input_term="abd-A", input_filenames=["../res/flymine_id_list_4.tsv"], 
                        preferences=Preferences(prefer_small_gene_names=True, prefer_first_selection_on_multiple=True, prefer_remember_selection=True, prefer_output=True, prefer_gene_db_id=True)))
    



# In[ ]:



