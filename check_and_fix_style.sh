echo "Black Fix"
black smiter tests
echo
echo 'Isort Fix'
isort -rc smiter tests
echo
# echo 'MyPy Check'
# mypy --ignore-missing-imports --pretty -p smiter
# echo
