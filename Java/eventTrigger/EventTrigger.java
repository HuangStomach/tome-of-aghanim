@Target({ElementType.Type, ElementType.METHOD})
@Retention(RetentionPolicy.RUNTIME)
public @interface EventTrigger {
    String value() default "";
    int order() default 100;
}

