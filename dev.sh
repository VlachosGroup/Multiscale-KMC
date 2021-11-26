#! /bin/bash

case $1 in

  "format")
    echo "Formatting code"
    find src include -type f \( -name "*.cpp" -o -name "*.hpp" \) -exec clang-format -i {} \;
  ;;
esac
