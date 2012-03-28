def get_match_dir(run, mock_catalog):
    d=deswl.files.collated_dir(run)
    d=os.path.join(d,'match-%s' % mock_catalog)
    if not os.path.exists(d):
        os.makedirs(d)
    return d

def get_match_files(run,mock_catalog):
    d=get_match_dir(run, mock_catalog)
    runf=os.path.join(d,'%s-radec.dat' % run)
    mockf=os.path.join(d,'%s-radec.dat' % mock_catalog)
    matchf='matches-%s-%s.dat' % (run,mock_catalog)
    matchf=os.path.join(d,matchf)
    scriptf=os.path.join(d,'domatch.sh')

    return {'dir':d,
            'runf':runf,
            'mockf':mockf,
            'matchf':matchf,
            'scriptf':scriptf}


