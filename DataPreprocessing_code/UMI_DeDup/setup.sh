
#try simly determine the pwd and write to the runMerge.sh script for correct directory path.

#determine the path
path=`pwd`

#determine the c++ library is available
echo "**check the neccessary libary.............."
FILE=./packages/alignutil_dir/libAlignUtil.so.1.0
if [ ! -f "$FILE" ]; then
    echo "$FILE does not exist!! Please compile to build it......"
	echo "setup quits............."
fi
echo "** Done.........."
echo "**Start writing the runMerge script........"
#echo 
echo "#!$(which bash)" > runMerge.sh
echo "#getting the correct directory for calling the python code" >>runMerge.sh
echo "${path}/merge_umi3_feng.py" >>runMerge.sh

echo "setup finished successfully..........."


