def filter_out_repeats(data_frame: pd.DataFrame, col_name: str) -> pd.DataFrame:
    """Filter out siRNA candidates with motifs of low sequence complexity such as “GGGG,” “UUUU,” “CCCC,” and “AAAA” """

    pattern_repeats = r'.*(A{5,16}|C{5,16}|G{5,16}|U{5,16}).*'
    return data_frame[~data_frame[col_name].str.contains(pattern_repeats)]


def filter_out_motifs(data_frame: pd.DataFrame, col_name: str) -> pd.DataFrame:
    """Filter out siRNA candidates that are vulnerable to nonspecific silencing efficiency due to the presence of immune-stimulating
       motifs. Three sequence motifs, “UGUGU,” “GUCCUUCAA,” and “AUCGAU(N)nGGGG,” should be included in an immune-motif avoidance
       list. In addition, U-rich sequences or sequences with biased nucleotide content, such as (G + U) >> (C + A) or with A + U or
       G + U rich motifs, should also be avoided."""

    avoid_motifs = '.*UGUGU.*|.*GUCCUUCAA.*|.*AUCGAU[ACGU]+GGGG.*'
    return data_frame[~data_frame[col_name].str.contains(avoid_motifs)]


def filter_out_tm(data_frame: pd.DataFrame, col_name: str, tm_start: float, tm_end: float) -> pd.DataFrame:
    """Filter out entries with melting temperature of siRNA duplex lower than 'tm_start' and higher than 'tm_end' """
    return data_frame[(data_frame[col_name].between(tm_start, tm_end))]