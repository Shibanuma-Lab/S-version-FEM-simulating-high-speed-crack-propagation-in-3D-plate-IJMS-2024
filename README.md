# arrest

====

# Overview
What is this code?

## Description


## Requirement

## Install
sfem_linear配下にarrestディレクトリを作成する。
```
git clone git@gitlab.com:morita/sfem_linear.git
git clone git@github.com:FuruhashiFumihito/arrest.git
cd arrest
pip install -r requirements.txt
```

## Usage
arrestディレクトリからsrc/main.pyを実行する。
```
python3 src/main.py
```
コマンドは以下。
```
--h, --help  Show help message and text.
--test_start  Specify first test number. (default: 1)
--test_end  Specify last test number. (default: 1)
--step_start  Specify first step number. (default: 0)
--step_end  Specify first test number. (default: 300)
--delete  If True, delete ALL files in /inputfiles and /results. (default: False)
--particular  If True, this program runs in the specified step. (default: False)
--is_jonly  If True, calculate only J-integral. (default: False)
--is_meshonly If True, calculate only meshs. (default: False)
```

## Demo

## Licence

[MIT](https://github.com/tcnksm/tool/blob/master/LICENCE)

## Reference
[1]User Name, 'Paper Titile' Conference Name pp.xx 20XX

[tcnksm](https://github.com/tcnksm)

cf. [how to write readme](https://deeeet.com/writing/2014/07/31/readme/)
