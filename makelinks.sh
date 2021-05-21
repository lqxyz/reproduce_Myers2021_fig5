topdir="$1"  # e.g., ./data/sherwood
dirs='ULI_MEDIUM_SAMPLE'
for dir in $dirs
do
    echo $dir
    if cd $topdir/$dir 2>/dev/null
    then
      file=$(ls *[0-9].job* | grep -v mean | grep -v last | tail -1)
      echo file=$file
      if [ "$file" ]
      then
        rm -rf $dir.lastpartial.joblib
        ln -s $file $dir.lastpartial.joblib
      fi
      file=$(ls *[0-9].job* | grep mean | grep -v last | tail -1)
      echo file=$file
      if [ "$file" ]
      then
        rm -rf $dir.lastmean.joblib
        ln -s $file $dir.lastmean.joblib
      fi
      if [ -s $dir.joblib ]
      then
        rm -rf $dir.lastmean.joblib
        ln -s $dir.joblib  $dir.lastmean.joblib
      fi
    fi
done

