.PHONY: clean All

All:
	@echo "----------Building project:[ 2D-SPH-Formulation - Debug ]----------"
	@"$(MAKE)" -f  "2D-SPH-Formulation.mk"
clean:
	@echo "----------Cleaning project:[ 2D-SPH-Formulation - Debug ]----------"
	@"$(MAKE)" -f  "2D-SPH-Formulation.mk" clean
