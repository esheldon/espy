rsync -av \
      --exclude "*svn*" \
      --exclude "*git*" \
      --exclude "*swp" \
      --exclude "*~" \
      --exclude "*pyc" \
      . $DES_MODULES_ROOT/install-config/
