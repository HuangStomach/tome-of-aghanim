import javax.servlet.*;

import Me;

@Component
public class MeFilter implements Filter {
	@Autowired
	private Me me;

	@Override
	public void doFilter(ServletRequest request, ServletResponse response, FilterChain chain) 
			throws IOException, ServletException {
		HttpServletRequest httpRequest = (HttpServletRequest) request;
		Strong token = httpRequest.getHeader("SomeHeader");

		if (token != null) {
			me.setUserId(token);
		}
		finally {
			me.destory();
		}
	}
}