
package com.example.multidatasource.config;

import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.boot.context.properties.ConfigurationProperties;
import org.springframework.boot.jdbc.DataSourceBuilder;
import org.springframework.boot.orm.jpa.EntityManagerFactoryBuilder;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.context.annotation.Primary;
import org.springframework.data.jpa.repository.config.EnableJpaRepositories;
import org.springframework.orm.jpa.JpaTransactionManager;
import org.springframework.orm.jpa.LocalContainerEntityManagerFactoryBean;
import org.springframework.transaction.PlatformTransactionManager;
import org.springframework.transaction.annotation.EnableTransactionManagement;

import javax.persistence.EntityManagerFactory;
import javax.sql.DataSource;
import java.util.HashMap;
import java.util.Map;

@Configuration
@EnableJpaRepositories(
    basePackages = "com.example.multidatasource.repository.c",
    entityManagerFactoryRef = "entityManagerFactoryC",
    transactionManagerRef = "transactionManagerC"
)
public class DataSourceConfig {
    @Bean(name = "dataSourceC")
    @ConfigurationProperties(prefix = "spring.datasource.c")
    public DataSource dataSourceC() {
        return DataSourceBuilder.create().build();
    }

    @Bean(name = "entityManagerFactoryC")
    public LocalContainerEntityManagerFactoryBean entityManagerFactoryC(
            EntityManagerFactoryBuilder builder,
            @Qualifier("dataSourceC") DataSource dataSource) {
        return builder
                .dataSource(dataSource)
                .packages("com.example.multidatasource.model.c")
                .persistenceUnit("c")
                .properties(jpaPropertiesC())
                .build();
    }

    private Map<String, Object> jpaPropertiesC() {
        Map<String, Object> props = new HashMap<>();
        props.put("hibernate.hbm2ddl.auto", "none");
        props.put("hibernate.dialect", "org.hibernate.dialect.MySQL8Dialect");
        return props;
    }

    @Bean(name = "transactionManagerC")
    public PlatformTransactionManager transactionManagerC(
            @Qualifier("entityManagerFactoryC") EntityManagerFactory entityManagerFactory) {
        return new JpaTransactionManager(entityManagerFactory);
    }
}
