def step2str(step: int) -> str:
    """
        整数を5桁かつ先頭0の文字列に変換する関数
        特にstepをfile名に埋め込む際に使う
        ex. step = 12 -> return "00012"
    """
    if 0 <= step <= 9:
        return "0000" + str(step)
    elif 10 <= step <= 99:
        return "000" + str(step)
    elif 100 <= step <= 999:
        return "00" + str(step)
    elif 1000 <= step <= 9999:
        return "0" + str(step)
    elif 10000 <= step <= 99999:
        return str(step)