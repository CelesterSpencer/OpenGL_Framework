//============================================================================
// Distributed under the MIT License. Author: Adrian Derstroff
//============================================================================

#version 430

in vec2 uv;
flat in int   id;

out vec3 FragColor;
layout (depth_less) out float gl_FragDepth; // Makes optimizations possible

void main()
{
    /*
     * render only sphere on the triangle billboard
     */
    float distance = length(uv);
    if(distance > 1.0)
    {
        discard;
    }

    // Output color
    FragColor = vec3(float(id), float(id), float(gl_PrimitiveID + 1));
}
