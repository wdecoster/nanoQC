cd ../
python -m build
cd ../../
mkdir testing-build      
cd testing-build      
python3 -m venv freshbuild               
. ./freshbuild/bin/activate          
cp -r ../nanoQC/dist .         
pip install dist/nanoQC-0.9.4-py3-none-any.whl

nanoQC -h 
cd ../
rm -r testing-build      