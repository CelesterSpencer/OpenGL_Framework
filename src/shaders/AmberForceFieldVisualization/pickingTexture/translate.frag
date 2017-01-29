//============================================================================
// Distributed under the MIT License. Author: Adrian Derstroff
//============================================================================

#version 430

in int axesID;

out vec3 FragColor;
layout (depth_less) out float gl_FragDepth; // Makes optimizations possible

void main() {
	// Output color
    FragColor = vec3(float(axesID), float(axesID), float(gl_PrimitiveID + 1));
}
