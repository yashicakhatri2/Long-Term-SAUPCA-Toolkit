function OUT = getJOrder1(G,H,J2,L,RE,b,cs,g,h,hs,ks,mu,sSun)
%getJOrder1
%    OUT = getJOrder1(G,H,J2,L,RE,B,CS,g,h,HS,KS,MU,sSun)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    25-Apr-2022 10:08:27

t2 = cos(g);
t3 = cos(ks);
t4 = sin(g);
t5 = sin(ks);
t6 = G.^2;
t7 = G.^3;
t9 = H.^2;
t10 = H.^3;
t12 = J2.^2;
t13 = L.^2;
t14 = RE.^2;
t16 = cs+1.0;
t17 = g.*2.0;
t18 = mu.^4;
t19 = mu.^6;
t22 = 1.0./G;
t29 = 1.0./L;
t35 = cs-1.0;
t36 = -h;
t37 = -hs;
t38 = -ks;
t39 = 1.0./mu;
t8 = t6.^2;
t11 = t9.^2;
t15 = t14.^2;
t20 = cos(t17);
t21 = sin(t17);
t23 = 1.0./t6;
t24 = 1.0./t7;
t26 = t22.^5;
t28 = t22.^7;
t30 = 1.0./t13;
t31 = t29.^3;
t33 = t29.^5;
t43 = H.*t22;
t46 = g+h+ks+t37;
t47 = g+hs+ks+t36;
t76 = g+h+t37+t38;
t77 = g+hs+t36+t38;
t25 = 1.0./t8;
t27 = t23.^3;
t32 = t30.^2;
t34 = t30.^3;
t41 = t24.^3;
t42 = t23.^5;
t44 = H.*t23.*1.6e+1;
t45 = H.*t24.*3.2e+1;
t48 = t43+1.0;
t49 = t9.*t23;
t50 = t6.*t30;
t52 = H.*t23.*8.4e+1;
t53 = H.*t24.*1.68e+2;
t54 = t43-1.0;
t57 = cos(t46);
t58 = cos(t47);
t59 = sin(t46);
t60 = sin(t47);
t64 = t9.*t24.*1.6e+1;
t66 = t10.*t26.*4.0e+1;
t71 = t9.*t24.*8.4e+1;
t72 = t10.*t26.*1.68e+2;
t74 = t10.*t26.*2.88e+2;
t82 = cos(t76);
t83 = cos(t77);
t84 = sin(t76);
t85 = sin(t77);
t90 = t10.*t23.*t31.*8.64e+2;
t40 = t25.^2;
t51 = -t45;
t55 = t49.*3.0;
t56 = t50.*6.0;
t61 = -t49;
t62 = t49.*8.0;
t63 = t49.*1.5e+1;
t65 = t49.*4.2e+1;
t67 = -t50;
t68 = t49-1.0;
t69 = t50-1.0;
t73 = t9.*t25.*2.88e+2;
t75 = t11.*t27.*2.88e+2;
t106 = (t16.*t22.*t58)./4.0;
t107 = (t16.*t23.*t58)./4.0;
t108 = (t16.*t22.*t60)./4.0;
t111 = (t22.*t35.*t57)./4.0;
t112 = (t23.*t35.*t57)./4.0;
t113 = (t22.*t35.*t59)./4.0;
t116 = H.*J2.*t14.*t18.*t26.*t32.*(9.0./2.0);
t117 = H.*J2.*t14.*t18.*t27.*t31.*(1.5e+1./2.0);
t119 = (H.*t16.*t23.*t60)./4.0;
t129 = (H.*t23.*t35.*t59)./4.0;
t134 = J2.*t9.*t14.*t18.*t27.*t32.*(9.0./2.0);
t141 = (t16.*t22.*t82)./4.0;
t142 = (t16.*t23.*t82)./4.0;
t143 = (t16.*t22.*t84)./4.0;
t149 = (H.*t16.*t23.*t84)./4.0;
t153 = (t22.*t35.*t83)./4.0;
t155 = (t23.*t35.*t83)./4.0;
t157 = (t22.*t35.*t85)./4.0;
t159 = H.*t16.*t23.*t82.*(-1.0./4.0);
t161 = (H.*t23.*t35.*t85)./4.0;
t162 = (t35.*t48.*t57)./4.0;
t163 = (t16.*t54.*t58)./4.0;
t164 = (t35.*t48.*t59)./4.0;
t165 = (t16.*t54.*t60)./4.0;
t169 = H.*t23.*t35.*t83.*(-1.0./4.0);
t175 = (t16.*t48.*t84)./4.0;
t179 = (t16.*t48.*t82)./4.0;
t187 = (t35.*t54.*t85)./4.0;
t190 = (t35.*t54.*t83)./4.0;
t70 = -t62;
t78 = t61+1.0;
t79 = t55-1.0;
t80 = t67+1.0;
t81 = t56-6.0;
t86 = t63-1.0;
t87 = t68.^2;
t100 = H.*t23.*t68.*2.0e+1;
t101 = H.*t24.*t68.*4.0e+1;
t104 = t52.*t68;
t105 = t53.*t68;
t114 = -t108;
t115 = t9.*t24.*t68.*2.0e+1;
t118 = H.*t107;
t121 = -t113;
t126 = t68.*t71;
t127 = H.*t112;
t130 = -t119;
t140 = -t129;
t147 = -t134;
t148 = H.*t142;
t152 = -t141;
t154 = -t142;
t156 = -t143;
t160 = H.*t155;
t166 = G.*H.*t20.*t32.*t68.*3.6e+2;
t167 = -t153;
t168 = -t155;
t170 = -t161;
t171 = -t162;
t172 = -t163;
t173 = -t164;
t176 = G.*t9.*t20.*t33.*t68.*3.6e+2;
t177 = H.*t6.*t20.*t33.*t68.*3.6e+2;
t180 = -t175;
t186 = -t179;
t191 = -t187;
t199 = -t190;
t88 = t79.^2;
t89 = t87.*5.0;
t91 = sqrt(t78);
t92 = sqrt(t80);
t94 = t87.*2.1e+1;
t99 = t23.*t79.*4.8e+1;
t102 = -t100;
t109 = H.*t24.*t79.*9.6e+1;
t110 = H.*t31.*t79.*1.44e+2;
t123 = G.*H.*t32.*t79.*4.32e+2;
t124 = -t115;
t125 = t9.*t24.*t79.*4.8e+1;
t137 = t9.*t25.*t79.*1.44e+2;
t138 = t9.*t31.*t79.*1.44e+2;
t139 = t9.*t32.*t79.*4.32e+2;
t174 = G.*H.*t20.*t32.*t86.*2.4e+1;
t181 = G.*t9.*t20.*t33.*t86.*2.4e+1;
t182 = H.*t6.*t20.*t33.*t86.*2.4e+1;
t192 = t52+t104;
t200 = t10.*t20.*t24.*t30.*t81.*1.2e+2;
t202 = J2.*t14.*t18.*t25.*t32.*t79.*(9.0./4.0);
t207 = t71+t126;
t208 = H.*t20.*t30.*t68.*t81.*3.0e+1;
t209 = H.*t20.*t31.*t68.*t81.*6.0e+1;
t212 = H.*t20.*t30.*t81.*t86.*2.0;
t213 = H.*t20.*t31.*t81.*t86.*4.0;
t221 = t7.*t20.*t32.*t68.*t86.*1.2e+1;
t222 = t8.*t20.*t33.*t68.*t86.*1.2e+1;
t225 = t51+t66+t101;
t230 = t7.*t20.*t33.*t68.*t86.*7.2e+1;
t234 = t9.*t20.*t22.*t30.*t68.*t81.*3.0e+1;
t235 = t9.*t20.*t22.*t31.*t68.*t81.*6.0e+1;
t236 = t53+t72+t105;
t239 = t9.*t20.*t22.*t30.*t81.*t86.*2.0;
t240 = t9.*t20.*t22.*t31.*t81.*t86.*4.0;
t246 = G.*t20.*t30.*t68.*t81.*t86.*2.0;
t247 = G.*t20.*t31.*t68.*t81.*t86.*4.0;
t248 = t20.*t50.*t68.*t81.*t86;
t252 = t6.*t20.*t31.*t68.*t81.*t86.*2.0;
t254 = t20.*t67.*t68.*t81.*t86;
t268 = t113+t114+t156+t157;
t269 = t119+t140+t149+t170;
t272 = t162+t172+t186+t190;
t273 = t162+t172+t179+t199;
t274 = t165+t173+t175+t191;
t93 = t88.*4.0;
t95 = 1.0./t91;
t97 = 1.0./t92;
t103 = H.*t99;
t120 = -t110;
t122 = G.*t110;
t128 = sSun.*t2.*t3.*t91;
t131 = sSun.*t2.*t5.*t91;
t132 = sSun.*t3.*t4.*t91;
t133 = sSun.*t4.*t5.*t91;
t135 = -t123;
t136 = -t125;
t144 = t7.*t31.*t88.*1.2e+1;
t145 = t6.*t31.*t88.*3.6e+1;
t146 = t7.*t32.*t88.*3.6e+1;
t150 = -t138;
t151 = -t139;
t158 = t6.*t32.*t88.*1.08e+2;
t178 = t70+t89;
t188 = t44+t102;
t198 = t65+t94-2.2e+1;
t203 = -t200;
t204 = -t202;
t205 = t64+t124;
t206 = G.*t30.*(t62-t89).*-6.0;
t211 = t50.*(t62-t89).*-3.0;
t214 = -t208;
t215 = t7.*t33.*(t62-t89).*-1.2e+1;
t216 = -t212;
t219 = G.*t30.*t192.*8.0;
t227 = t50.*t192.*4.0;
t228 = -t221;
t229 = -t222;
t232 = t6.*t31.*t192.*8.0;
t233 = -t230;
t241 = t50.*t207.*4.0;
t243 = t6.*t31.*t207.*8.0;
t245 = t6.*t31.*t69.*(t62-t89).*-6.0;
t249 = t50.*t225.*3.0;
t250 = -t246;
t251 = -t247;
t255 = -t252;
t256 = t50.*t236.*4.0;
t278 = b.*t13.*t39.*t92.*(t108+t121+t143-t157).*(-3.0./2.0);
t281 = b.*t13.*t39.*t92.*t269.*(3.0./2.0);
t290 = L.*b.*t39.*t92.*(t164-t165+t180+t187).*-3.0;
t292 = b.*t13.*t39.*t92.*t272.*(3.0./2.0);
t293 = b.*t13.*t39.*t92.*t273.*(3.0./2.0);
t96 = t95.^3;
t98 = t97.^3;
t183 = H.*sSun.*t2.*t5.*t23.*t95;
t184 = H.*sSun.*t3.*t4.*t23.*t95;
t185 = H.*sSun.*t4.*t5.*t23.*t95;
t189 = H.*sSun.*t4.*t5.*t24.*t95.*2.0;
t193 = sSun.*t2.*t5.*t9.*t24.*t95;
t194 = sSun.*t3.*t4.*t9.*t24.*t95;
t195 = sSun.*t4.*t5.*t9.*t24.*t95;
t210 = G.*t30.*t188.*6.0;
t217 = G.*t30.*t198.*8.0;
t218 = G.*t31.*t198.*1.6e+1;
t220 = t50.*t188.*3.0;
t223 = -t219;
t224 = t50.*t198.*4.0;
t231 = t6.*t31.*t198.*8.0;
t237 = -t232;
t238 = t50.*t205.*3.0;
t242 = -t241;
t244 = -t243;
t253 = t6.*t31.*t69.*t188.*6.0;
t257 = t93+t211;
t285 = t128+t162+t163+t179+t190;
t287 = t132+t164+t165+t175+t187;
t288 = t133+t171+t172+t179+t190;
t289 = t131+t164+t165+t180+t191;
t291 = G.*b.*t39.*t97.*(t164-t165+t180+t187).*(-3.0./2.0);
t294 = b.*t6.*t29.*t39.*t97.*(t164-t165+t180+t187).*(-3.0./2.0);
t196 = sSun.*t4.*t5.*t10.*t26.*t96;
t197 = -t184;
t201 = -t194;
t226 = -t220;
t258 = G.*t30.*t257.*2.0;
t259 = G.*t31.*t257.*4.0;
t261 = t6.*t31.*t257.*2.0;
t262 = t69.*t257;
t263 = G.*t30.*(t220-H.*t23.*t79.*4.8e+1).*-2.0;
t264 = G.*t30.*(t220-H.*t23.*t79.*4.8e+1).*2.0;
t265 = t6.*t31.*(t220-H.*t23.*t79.*4.8e+1).*-2.0;
t266 = t6.*t31.*(t220-H.*t23.*t79.*4.8e+1).*2.0;
t267 = -t69.*(t220-H.*t23.*t79.*4.8e+1);
t270 = t136+t206+t238;
t271 = t6.*t31.*(t125-t238+G.*t30.*(t62-t89).*6.0).*-2.0;
t275 = -t69.*(t125-t238+G.*t30.*(t62-t89).*6.0);
t276 = t74+t109+t210+t249;
t277 = t106+t111+t152+t167+t185;
t279 = t114+t121+t143+t157+t183;
t282 = t118+t127+t159+t169+t195;
t283 = t130+t140+t149+t161+t193;
t306 = L.*b.*t39.*t92.*t287.*3.0;
t307 = L.*b.*t39.*t92.*t289.*3.0;
t308 = G.*b.*t39.*t97.*t287.*(3.0./2.0);
t309 = b.*t13.*t39.*t92.*t285.*(3.0./2.0);
t311 = G.*b.*t39.*t97.*t289.*(3.0./2.0);
t313 = b.*t13.*t39.*t92.*t288.*(3.0./2.0);
t315 = b.*t6.*t29.*t39.*t97.*t287.*(3.0./2.0);
t316 = b.*t7.*t31.*t39.*t98.*t288.*(3.0./2.0);
t317 = b.*t6.*t29.*t39.*t97.*t289.*(3.0./2.0);
t260 = t103+t226;
t280 = t108+t113+t143+t157+t197;
t284 = t119+t129+t149+t161+t201;
t286 = t69.*t276;
t295 = L.*b.*t39.*t92.*t277.*3.0;
t296 = G.*b.*t39.*t97.*t277.*(3.0./2.0);
t297 = b.*t13.*t39.*t92.*t279.*(3.0./2.0);
t300 = L.*b.*t39.*t92.*t282.*3.0;
t302 = t107+t112+t154+t168+t189+t196;
t303 = b.*t6.*t29.*t39.*t97.*t277.*(3.0./2.0);
t304 = b.*t13.*t39.*t92.*t283.*(3.0./2.0);
t310 = b.*t6.*t29.*t39.*t97.*t282.*(3.0./2.0);
t314 = -t313;
t318 = -t316;
t321 = t93+t144+t224+t254+t262;
t324 = t103+t122+t214+t216+t227+t267;
t327 = t146+t229+t231+t245+t255+t261;
t330 = t135+t177+t182+t209+t213+t237+t253+t266;
t335 = t136+t145+t150+t217+t228+t234+t239+t242+t250+t258+t275;
t336 = t12.*t15.*t19.*t30.*t41.*(t125+t138-t145-t217+t221-t234-t239+t241+t246-t258+t69.*(t125-t238+G.*t30.*(t62-t89).*6.0)).*(-1.0./1.28e+2);
t298 = b.*t13.*t39.*t92.*t280.*(3.0./2.0);
t299 = -t296;
t301 = -t300;
t305 = b.*t13.*t39.*t92.*t284.*(3.0./2.0);
t312 = -t310;
t319 = b.*t13.*t39.*t92.*t302.*(3.0./2.0);
t322 = t12.*t15.*t19.*t30.*t42.*t321.*(9.0./1.28e+2);
t325 = (t12.*t15.*t19.*t30.*t41.*t324)./1.28e+2;
t326 = t12.*t15.*t19.*t29.*t42.*t324.*(9.0./1.28e+2);
t328 = t12.*t15.*t19.*t29.*t42.*t327.*(9.0./1.28e+2);
t331 = (t12.*t15.*t19.*t29.*t41.*t330)./1.28e+2;
t333 = t74+t90+t109+t120+t166+t174+t203+t223+t256+t264+t286;
t320 = -t319;
t323 = -t322;
t329 = -t328;
t332 = -t331;
t334 = (t12.*t15.*t19.*t29.*t41.*t333)./1.28e+2;
t337 = t116+t295+t303+t325+t332;
t338 = t117+t299+t320+t326+t334;
et1 = b.*t39.*t97.*t288.*(3.0./2.0)+G.*b.*t39.*t97.*t282.*3.0+b.*t39.*t50.*t98.*t288.*(3.0./2.0)+b.*t13.*t39.*t92.*((H.*t16.*t24.*t58)./2.0+(H.*t24.*t35.*t57)./2.0-(H.*t16.*t24.*t82)./2.0-(H.*t24.*t35.*t83)./2.0+sSun.*t4.*t5.*t9.*t25.*t95.*3.0+sSun.*t4.*t5.*t11.*t27.*t96).*(3.0./2.0)-t12.*t15.*t19.*t22.^11.*t29.*t321.*(4.5e+1./6.4e+1)-t12.*t15.*t19.*t29.*t42.*(t125+t138-t145-t217+t221-t234-t239+t241+t246-t258+t69.*(t125-t238+G.*t30.*(t62-t89).*6.0)).*(9.0./6.4e+1);
et2 = t12.*t15.*t19.*t29.*t41.*(t75+t137+t30.*t198.*8.0+t30.*t257.*2.0+t50.*(t9.*t25.*2.52e+2+t11.*t27.*1.68e+2+t9.*t25.*t68.*2.52e+2).*4.0+t69.*(t75+t137-t30.*(t62-t89).*6.0+t50.*(t9.*t25.*-4.8e+1+t11.*t27.*4.0e+1+t9.*t25.*t68.*6.0e+1).*3.0+G.*t30.*t205.*1.2e+1)-G.*t30.*(t125-t238+G.*t30.*(t62-t89).*6.0).*4.0+G.*t31.*t88.*7.2e+1-G.*t30.*t207.*1.6e+1+t11.*t24.*t31.*8.64e+2+t9.*t20.*t32.*t68.*7.2e+2-t9.*t22.*t31.*t79.*4.32e+2+t9.*t20.*t32.*t86.*4.8e+1-t11.*t20.*t25.*t30.*t81.*1.2e+2-t6.*t20.*t32.*t68.*t86.*6.0e+1+t20.*t30.*t49.*t68.*t81.*3.0e+1+t20.*t30.*t49.*t81.*t86.*2.0-t20.*t30.*t68.*t81.*t86.*2.0).*(-1.0./1.28e+2);
et3 = J2.*t9.*t14.*t18.*t28.*t31.*(-2.7e+1./2.0)-J2.*t14.*t18.*t26.*t31.*t79.*3.0;
mt1 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t307-t317+(t12.*t15.*t19.*t29.*t41.*(t8.*t21.*t33.*t68.*t86.*2.4e+1+t6.*t21.*t31.*t68.*t81.*t86.*4.0))./1.28e+2+(t12.*t15.*t19.*t21.*t28.*t32.*t68.*t81.*t86)./6.4e+1,-t304+t311-(t12.*t15.*t19.*t29.*t41.*(t7.*t21.*t32.*t68.*t86.*2.4e+1-t9.*t21.*t22.*t30.*t68.*t81.*6.0e+1-t9.*t21.*t22.*t30.*t81.*t86.*4.0+G.*t21.*t30.*t68.*t81.*t86.*4.0))./1.28e+2+t12.*t15.*t19.*t21.*t31.*t40.*t68.*t81.*t86.*(9.0./6.4e+1),t297-(t12.*t15.*t19.*t29.*t41.*(H.*t21.*t30.*t68.*t81.*6.0e+1+H.*t21.*t30.*t81.*t86.*4.0))./1.28e+2,0.0,t314+(t12.*t15.*t19.*t20.*t28.*t31.*t68.*t81.*t86)./3.2e+1,t292,0.0,t309,t290+t294];
mt2 = [-t281+G.*b.*t39.*t97.*(t164-t165+t180+t187).*(3.0./2.0),b.*t13.*t39.*t92.*(t108+t121+t143-t157).*(3.0./2.0),0.0,t292,b.*t13.*t39.*t92.*(t162+t163+t186+t199).*(3.0./2.0),0.0,t293,t32.*1.0./t39.^2.*-3.0-b.*t39.*t92.*t288.*3.0-b.*t39.*t50.*t97.*t288.*(3.0./2.0)+b.*t8.*t32.*t39.*t98.*t288.*(3.0./2.0)-(t12.*t15.*t19.*t31.*t41.*t321)./6.4e+1-(t12.*t15.*t19.*t30.*t41.*t327)./6.4e+1+(t12.*t15.*t19.*t29.*t41.*(t7.*t33.*t88.*-1.44e+2-t6.*t32.*t198.*2.4e+1-t6.*t32.*t257.*6.0+t8.*t34.*(t62-t89).*2.4e+1+t6.*t32.*t69.*(t62-t89).*1.8e+1+t8.*t20.*t34.*t68.*t86.*8.4e+1+t6.*t20.*t32.*t68.*t81.*t86.*6.0))./1.28e+2-J2.*t14.*t18.*t24.*t33.*t79.*3.0];
mt3 = [t147+t204+t301+t312+t318+t323+t329+t336+(t12.*t15.*t19.*t29.*t41.*(t151+t158+t176+t181+t215+t218+t233+t235+t240+t244+t251+t259+t271+t69.*(t6.*t31.*t205.*6.0-G.*t31.*(t62-t89).*1.2e+1)))./1.28e+2,t337,0.0,t307+t317-t12.*t15.*t19.*t21.*t26.*t34.*t68.*t86.*(3.0./1.6e+1)-t12.*t15.*t19.*t21.*t28.*t32.*t68.*t81.*t86.*(3.0./6.4e+1),L.*b.*t39.*t92.*(t164-t165+t180+t187).*3.0+b.*t6.*t29.*t39.*t97.*(t164-t165+t180+t187).*(3.0./2.0),0.0,t306+t315,t147+t204+t301+t312+t318+t323+t329+t336+(t12.*t15.*t19.*t29.*t41.*(t151+t158+t176+t181+t215+t218+t233+t235+t240+t244+t251+t259+t271+t6.*t31.*t69.*t205.*6.0-G.*t31.*t69.*(t62-t89).*1.2e+1))./1.28e+2,et1+et2+et3,t338,0.0];
mt4 = [t304-t311+t12.*t15.*t19.*t21.*t27.*t33.*t68.*t86.*(3.0./1.6e+1)-t9.*t12.*t15.*t19.*t21.*t31.*t42.*t68.*t81.*(1.5e+1./3.2e+1)-(t9.*t12.*t15.*t19.*t21.*t31.*t42.*t81.*t86)./3.2e+1-t12.*t15.*t19.*t21.*t31.*t40.*t68.*t81.*t86.*(7.0./6.4e+1),t281+t291,0.0,-t305-t308,t337,t338];
mt5 = [b.*t13.*t39.*t92.*(sSun.*t4.*t5.*t23.*t95+sSun.*t4.*t5.*t9.*t25.*t96).*(3.0./2.0)-J2.*t14.*t18.*t26.*t31.*(3.0./2.0)-(t12.*t15.*t19.*t29.*t41.*(t73+t99+t50.*(t23.*8.4e+1+t9.*t25.*1.68e+2+t23.*t68.*8.4e+1).*4.0+t69.*(t73+t99+t50.*(t23.*-1.6e+1+t9.*t25.*4.0e+1+t23.*t68.*2.0e+1).*3.0)+G.*t31.*t79.*1.44e+2+t9.*t22.*t31.*8.64e+2-t20.*t30.*t49.*t81.*1.2e+2-t20.*t30.*t68.*t81.*3.0e+1-t20.*t30.*t81.*t86.*2.0))./1.28e+2,0.0,-t297+H.*t12.*t15.*t19.*t21.*t31.*t41.*t68.*t81.*(1.5e+1./3.2e+1)+(H.*t12.*t15.*t19.*t21.*t31.*t41.*t81.*t86)./3.2e+1,t278,0.0,t298,-t306-t315,t305+t308,-t298,0.0,t309,t293,0.0,t314,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
OUT = reshape([mt1,mt2,mt3,mt4,mt5],8,8);
