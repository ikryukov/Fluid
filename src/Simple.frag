#version 410 core
in vec4 fColor;
in vec3 fNormal;

vec3 Light = vec3(0.0, 4.0, 7.0);
vec4 Color = vec4(0.2, 0.4, 0.5, 1.0);

out vec4 fragColor;

void main(void)
{
	const float fAmbient = 0.4;
	Light = normalize(Light);
	//float depth = (fShadowMapCoord.z / fShadowMapCoord.w);
	//float depth_light = texture2DProj(shadowMapTex, fShadowMapCoord).r;
	//highp float visibility = depth <= depth_light ? 1.0 : 0.2;
	//float visibility = clamp(step(depth, depth_light) + 0.2, 0.0, 1.0);
	float visibility = 1.0;
	fragColor = fColor * max(0.0, dot(fNormal, Light)) * visibility;
	//gl_FragColor = fColor;
}
