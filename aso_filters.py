def filter_out_repeats(data_frame: pd.DataFrame, col_name: str) -> pd.DataFrame:
    """
    Filter out ASO candidates with motifs of low sequence complexity with 5-16 repeating nucleosides
    e.g. “GGGGG,” “TTTTT,” “CCCCC,” and “AAAAA”, "TGTGTG", "ACACAC" etc.
    Return dataframe containing entries without such motifs
    """
    pattern_repeats = r'.*(A{5,16}|C{5,16}|G{5,16}|T{5,16}|(AC){3,8}|(AC){3,8}|(AG){3,8}|(AT){3,8}|(CA){3,8}|(CG){3,8}|(CT){3,8}|(GA){3,8}|(GC){3,8}|(GT){3,8}|(TA){3,8}|(TC){3,8}|(TG){3,8}).*'
    return data_frame[~data_frame[col_name].str.contains(pattern_repeats)]


def filter_gc_content(data_frame: pd.DataFrame, col_name: str, min_gc: float, max_gc: float) -> pd.DataFrame:
    """
    Filter dataframe according to GC-content.
    Return dataframe containing entries with GC-content between 'min_gc' and 'max_gc'
    """
    return data_frame[data_frame[col_name].between(min_gc, max_gc)]


def filter_pLfold(data_frame: pd.DataFrame, col_name: str, pLfold_threshold: float) -> pd.DataFrame:
    """
    Filter dataframe according to pLfold value.
    Return dataframe containing entries with pLfold value higher or equal than 'pLfold_threshold'
    """
    return data_frame[data_frame[col_name] >= pLfold_threshold]


def filter_ASO_ASO_freeenergy(data_frame: pd.DataFrame, col_name: str, ASO_ASO_threshold: float) -> pd.DataFrame:
    """
    Filter dataframe according to ASO_ASO_freeenergy value.
    Return dataframe containing entries with ASO_ASO_freeenergy value higher or equal than 'ASO_ASO_threshold'
    """
    return data_frame[data_frame[col_name] >= ASO_ASO_threshold]


def filter_ASO_transcriptfreeenergy(data_frame: pd.DataFrame, col_name: str, min_tfe: float, max_tfe: float) -> pd.DataFrame:
    """
    Filter dataframe according to ASO Transcript free energy.
    Return dataframe containing entries with ASO Transcript free energy between 'min_tfe' and 'max_tfe'
    """
    return data_frame[data_frame[col_name].between(min_tfe, max_tfe)]